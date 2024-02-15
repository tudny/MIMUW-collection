module InscJvm where

import AbstractSyntaxTree
import Control.Monad.Reader
import Control.Monad.State
import Data.Map ((!))
import Exec (exec)
import JvmUtils
import System.Environment
import System.FilePath
import System.Process (system)
import Utils (decodeFilename)

type PM = StateT [Code] (ReaderT Variables IO)

main :: IO ()
main =
  exec
    go
    [ "Usage: insc_jvm [file.ins]",
      "Compiles the given instant program to JVM bytecode.",
      "The output is written to file.j and file.class."
    ]

go :: String -> AST -> IO ()
go source ast = do
  let (path, name, _) = decodeFilename source
  buildCode path name ast

buildCode :: String -> String -> AST -> IO ()
buildCode path classname ast = do
  let vars = getVariables ast
  code <- parseAST ast
  let stackSize' = mesureStack code
  let localsSize = getCount vars
  jasminFile <- pareseJvmTemplate classname stackSize' localsSize code
  saveJasminFile path classname jasminFile
  runJasmin path classname

runJasmin :: String -> String -> IO ()
runJasmin path classname = do
  execPath <- getExecutablePath
  let execBasePath = takeDirectory execPath
  let jasminPath = execBasePath ++ "/lib/jasmin.jar"
  let filename = makeJasminFilePath path classname
  let command = makeJasminCommand jasminPath filename path
  _ <- system command
  pure ()

makeJasminCommand :: String -> String -> String -> String
makeJasminCommand jasminPath filename "" = "java -jar " ++ jasminPath ++ " " ++ filename
makeJasminCommand jasminPath filename path = "java -jar " ++ jasminPath ++ " " ++ filename ++ " -d " ++ path

makeJasminFilePath :: String -> String -> String
makeJasminFilePath "" name = name ++ ".j"
makeJasminFilePath path name = path ++ "/" ++ name ++ ".j"

mesureStack :: [Code] -> Integer
mesureStack [] = 0
mesureStack cs = maximum . scanl1 (+) . map stackImpact $ cs

stackImpact :: Code -> Integer
stackImpact (IPUSH _) = 1
stackImpact IADD = -1
stackImpact ISUB = -1
stackImpact IMUL = -1
stackImpact IDIV = -1
stackImpact (ISTORE _) = -1
stackImpact (ILOAD _) = 1
stackImpact IRETURN = -1
stackImpact IPRINT = -1
stackImpact SWAP = 0

parseAST :: AST -> IO [Code]
parseAST ast = do
  let vars = getVariables ast
  let stmts = getStmts ast
  (_, code) <- runReaderT (runStateT (codifyStmts stmts) []) vars
  pure $ reverse code

codifyStmts :: [IStmt] -> PM ()
codifyStmts stmts = do
  mapM_ codifyStmt stmts

codifyStmt :: IStmt -> PM ()
codifyStmt (IExpr expr) = do
  codifyExpr expr
  appendCodeM IPRINT
codifyStmt (IAss name expr) = do
  codifyExpr expr
  vars <- ask
  let index = getVars vars ! name
  appendCodeM $ ISTORE index

codifyExpr :: EExpr -> PM ()
codifyExpr (EExpr' (EInt c) _ _) = do
  appendCodeM $ IPUSH c
codifyExpr (EExpr' (EVar name) _ _) = do
  vars <- ask
  let index = getVars vars ! name
  appendCodeM $ ILOAD index
codifyExpr (EExpr' (EAdd e1 e2) _ s) = do
  codifyTwoBySize s False e1 e2
  appendCodeM IADD
codifyExpr (EExpr' (ESub e1 e2) _ s) = do
  codifyTwoBySize s True e1 e2
  appendCodeM ISUB
codifyExpr (EExpr' (EMul e1 e2) _ s) = do
  codifyTwoBySize s False e1 e2
  appendCodeM IMUL
codifyExpr (EExpr' (EDiv e1 e2) _ s) = do
  codifyTwoBySize s True e1 e2
  appendCodeM IDIV

codifyTwoBySize :: Bool -> Bool -> EExpr -> EExpr -> PM ()
codifyTwoBySize swap appendSwap e1 e2 = do
  if swap
    then do
      codifyExpr e2
      codifyExpr e1
      when appendSwap $ appendCodeM SWAP
    else do
      codifyExpr e1
      codifyExpr e2

appendCodeM :: Code -> PM ()
appendCodeM = modify . (:)
