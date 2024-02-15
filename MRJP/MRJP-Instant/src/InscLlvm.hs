module InscLlvm where

import AbstractSyntaxTree
import Control.Monad.Reader
import Control.Monad.State
import Data.Bifunctor
import Data.Map (keys)
import Exec (exec)
import LlvmUtils
import System.Environment
import System.FilePath
import System.Process
import Utils

type PM = StateT (Integer, [Code]) (ReaderT Variables IO)

main :: IO ()
main =
  exec
    go
    [ "Usage: insc_llvm [file.ins]",
      "Compiles the given instant program to LLVM IR.",
      "The output is written to file.ll and file.bc."
    ]

go :: String -> AST -> IO ()
go source ast = do
  let (path, name, _) = decodeFilename source
  buildCode path name ast

buildCode :: String -> String -> AST -> IO ()
buildCode path name ast = do
  code <- parseAST ast
  jjvmFile <- parseLlvmTemplate code
  saveLlvmFile path name jjvmFile
  runLlvm path name

runLlvm :: String -> String -> IO ()
runLlvm path name = do
  execPath <- getExecutablePath
  let execBasePath = takeDirectory execPath
  let filename_wo_ext = joinTwoStrWith path name "/"
  let command = "llvm-as -o " ++ filename_wo_ext ++ ".bc " ++ filename_wo_ext ++ ".ll"
  _ <- system command
  let command' = "llvm-link -o " ++ filename_wo_ext ++ ".bc " ++ filename_wo_ext ++ ".bc " ++ execBasePath ++ "/lib/runtime.bc"
  _ <- system command'
  putStrLn $ "Generated " ++ filename_wo_ext ++ ".bc"

joinTwoStrWith :: String -> String -> String -> String
joinTwoStrWith "" s2 _ = s2
joinTwoStrWith s1 "" _ = s1
joinTwoStrWith s1 s2 s = s1 ++ s ++ s2

parseAST :: AST -> IO [Code]
parseAST ast = do
  let vars = getVariables ast
  let stmts = getStmts ast
  allocas <- generateAllocas vars
  (_, (_, code)) <- runReaderT (runStateT (codifyStmts stmts) (0, allocas)) vars
  pure $ reverse code

generateAllocas :: Variables -> IO [Code]
generateAllocas vars = do
  let vars' = keys $ getVars vars
  let allocas = map (IALLOCA . ("%var" ++)) vars'
  pure allocas

codifyStmts :: [IStmt] -> PM ()
codifyStmts stmts = do
  mapM_ codifyStmt stmts

codifyStmt :: IStmt -> PM ()
codifyStmt (IAss v e) = do
  val <- codifyExpr e
  let reg = "%var" ++ v
  appendCodeM $ ISTORE val reg
codifyStmt (IExpr e) = do
  val <- codifyExpr e
  appendCodeM $ IPRINT val

codifyExpr :: EExpr -> PM String
codifyExpr (EExpr' (EInt i) _ _) = pure $ show i
codifyExpr (EExpr' (EVar v) _ _) = do
  resultRegNum <- incStateM
  let resultReg = "%temp" ++ show resultRegNum
  appendCodeM $ ILOAD resultReg $ "%var" ++ v
  pure resultReg
codifyExpr (EExpr' (EAdd e1 e2) _ _) = codifyTwoArg e1 e2 IADD
codifyExpr (EExpr' (ESub e1 e2) _ _) = codifyTwoArg e1 e2 ISUB
codifyExpr (EExpr' (EMul e1 e2) _ _) = codifyTwoArg e1 e2 IMUL
codifyExpr (EExpr' (EDiv e1 e2) _ _) = codifyTwoArg e1 e2 IDIV

codifyTwoArg :: EExpr -> EExpr -> (String -> String -> String -> Code) -> PM String
codifyTwoArg e1 e2 op = do
  val1 <- codifyExpr e1
  val2 <- codifyExpr e2
  resultRegNum <- incStateM
  let resultReg = "%temp" ++ show resultRegNum
  appendCodeM $ op resultReg val1 val2
  pure resultReg

incStateM :: PM Integer
incStateM = do
  modify $ Data.Bifunctor.first (+ 1)
  gets fst

appendCodeM :: Code -> PM ()
appendCodeM code = do
  modify $ Data.Bifunctor.second (code :)
  pure ()
