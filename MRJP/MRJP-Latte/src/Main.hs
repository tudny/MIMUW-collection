module Main where

import Backend.Core
import Control.Monad
import Frontend.TypeChecker
import Frontend.TypeErrors
import Frontend.Commons
import Latte.Abs
import Latte.Par
import System.Environment
import System.Exit
import Prelude


handleError :: String -> IO ()
handleError err = do
  putStrLnStderr "ERROR"
  putStrLnStderr err

tabFixedSize :: Int
tabFixedSize = 4

fixTabs :: String -> String
fixTabs = concatMap removeTabs
  where
    removeTabs :: Char -> String
    removeTabs '\t' = replicate tabFixedSize ' '
    removeTabs x = [x]

runFile :: FilePath -> IO ()
runFile filename = do
  putStrLn $ "Compiling " ++ filename ++ "..."
  content <- readFile filename
  runLexerParser filename $ fixTabs content

runLexerParser :: FilePath -> String -> IO ()
runLexerParser inputFile sourceCode = do
  putStrLn "Parsing..."
  let ts = myLexer sourceCode
  case pProgram ts of
    Left err -> do
      handleError err
      exitFailure
    Right tree -> do
      putStrLn "Parse successful!"
      runTypeChecker inputFile sourceCode tree

runTypeChecker :: FilePath -> String -> Program -> IO ()
runTypeChecker inputFile sourceCode tree = do
  putStrLn "Type checking..."
  tcResult <- typeChecker tree
  case tcResult of
    Left err -> do
      handleError "Type checking failed..."
      printErrHolder sourceCode err
      exitFailure
    Right (ast, astEnv) -> do
      putStrLn "Type checking successful!"
      runBackend inputFile ast astEnv
      putStrLnStderr "OK"
      exitSuccess

usage :: IO ()
usage = do
  putStrLn $
    unlines
      [ "usage: Call with one of the following argument combinations:",
        "  --help         Display this help message.",
        "  <file>         Compile content of <file>."
      ]

main :: IO ()
main = do
  args <- getArgs
  case args of
    ["--help"] -> usage >> exitSuccess
    [file] -> runFile file
    _ -> usage >> exitFailure
