module Exec where

import AbstractSyntaxTree (AST, parseTree)
import Distribution.Fields.Lexer (Token)
import Instant.Par (myLexer, pProgram)
import System.Environment (getArgs)
import System.Exit (exitFailure)
import Utils (exitWithCode)

type Compiler = String -> AST -> IO ()

type Usage = [String]

type Err = Either String

type ParseFun a = [Token] -> Err a

exec :: Compiler -> Usage -> IO ()
exec backend usage = do
  args <- getArgs
  case args of
    ["--help"] -> showUsage usage 0
    [file] -> runFile file backend
    _ -> showUsage usage 1

showUsage :: Usage -> Int -> IO ()
showUsage message code = do
  putStrLn $ unlines message
  exitWithCode code

runFile :: String -> Compiler -> IO ()
runFile file backend = do
  readFile file >>= run backend file

run :: Compiler -> String -> String -> IO ()
run backend filename code =
  case pProgram ts of
    Left err -> do
      putStrLn "\n Parse Failed...\n"
      putStrLn err
      exitFailure
    Right tree -> do
      case parseTree tree of
        Left err -> do
          putStrLn "\n Static Analysis Failed...\n"
          putStrLn err
          exitFailure
        Right ast -> do
          backend filename ast
  where
    ts = myLexer code
