module Main where

import System.Environment ( getArgs )
import System.Exit        ( exitFailure, exitWith, ExitCode (ExitSuccess, ExitFailure) )

import Src.Jabba.Abs   (Program)
import Src.Jabba.Lex   ( Token )
import Src.Jabba.Par   ( pProgram, myLexer )
import Src.Jabba.Skel  ()
import Src.TypeChecker ( typeCheck )
import Src.Evaluator   ( evaluate )
import Src.Errors      ( ErrHolder (ParserErr, ControlledExit) )
import Src.Utils       ( left )

type Err        = Either String
type ParseFun a = [Token] -> Err a

runFile :: ParseFun Program -> FilePath -> IO ()
runFile p f = putStrLn f >> readFile f >>= run p

run :: ParseFun Program -> String -> IO ()
run p s = case go of
  Left err -> do
    putStrLn "\nPrerun Error!\n"
    print err
    exitFailure
  Right prog -> do
    e <- evaluate prog s
    case e of
      Left (ControlledExit c) -> do
        putStrLn $ "Program exited with code " ++ show c
        exitWith $ if c == 0
           then ExitSuccess
           else ExitFailure c
      Left eh -> do
        putStrLn "\nRuntime Error!\n"
        print eh
        exitFailure
      Right () -> pure ()
  where
    ts = myLexer s
    go = do
      tree <- left ParserErr $ p ts
      typeCheck tree
      pure tree

usage :: IO ()
usage = do
  putStrLn $ unlines
    [ "usage: Call with one of the following argument combinations:"
    , "  --help          Display this help message."
    , "  (no arguments)  Parse stdin."
    , "  (filename)         Parse content of a file."
    ]

main :: IO ()
main = do
  args <- getArgs
  case args of
    ["--help"] -> usage
    []         -> getContents >>= run pProgram
    [f]        -> runFile pProgram f
    _          -> usage

