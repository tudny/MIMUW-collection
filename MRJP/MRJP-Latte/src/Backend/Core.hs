module Backend.Core where

import Backend.Codegen (codegen)
import qualified Backend.LLVMCode as L
import Backend.Renaming
import Data.List (intercalate)
import Frontend.AbstractSyntaxTree
import GHC.IO.Exception
import System.Environment
import System.FilePath
import System.Process
import Text.Pretty.Simple (pPrint)
import Utils (decodeFilename, exitWithCode, joinTwoStrWith)

runBackend :: FilePath -> AST -> ASTEnv -> IO ()
runBackend inputFile ast' astEnv' = do
  (ast, astEnv) <- rename ast' astEnv'
  pPrint ast
  pPrint astEnv
  program <- codegen ast astEnv
  putStrLn "====================================================="
  pPrint program
  putStrLn "====================================================="
  prelude <- readPrelude
  let llvmCode = prelude ++ L.codegen program
  putStrLn llvmCode
  putStrLn "====================================================="
  putStrLn "Running backend..."
  writeLLVMCode inputFile llvmCode

runtimeBC :: String
runtimeBC = "/lib/runtime.bc"

preludeLL :: String
preludeLL = "/lib/prelude.ll"

readPrelude :: IO String
readPrelude = do
  execPath <- getExecutablePath
  let execBasePath = takeDirectory execPath
  let preludePath = execBasePath ++ preludeLL
  readFile preludePath

writeLLVMCode :: FilePath -> String -> IO ()
writeLLVMCode inputFile content = do
  let (path, name, _) = decodeFilename inputFile
  execPath <- getExecutablePath

  let execBasePath = takeDirectory execPath
  let filenameWoExt = joinTwoStrWith path name "/"
  let filenameWithLL = filenameWoExt ++ ".ll"
  let filenameWithBC = filenameWoExt ++ ".bc"

  writeFile filenameWithLL content
  let command = "llvm-as -o " ++ filenameWithBC ++ " " ++ filenameWithLL
  putStrLn $ "Running command: " ++ command
  llvmCompileCode <- system command
  case llvmCompileCode of
    ExitSuccess -> return ()
    ExitFailure n -> do
      putStrLn $ "LLVM compilation failed with code: " ++ show n
      exitWithCode n

  let command' = "llvm-link -o " ++ filenameWithBC ++ " " ++ filenameWithBC ++ " " ++ execBasePath ++ runtimeBC
  putStrLn $ "Running command: " ++ command'
  llvmLinkCode <- system command'
  case llvmLinkCode of
    ExitSuccess -> return ()
    ExitFailure n -> do
      putStrLn $ "LLVM linking failed with code: " ++ show n
      exitWithCode n

  putStrLn $ "Saved LLVM code to file: " ++ filenameWithBC

saveCodeToFile :: String -> String -> IO ()
saveCodeToFile fileName code = do
  writeFile fileName code
  putStrLn $ "Saved code to file: " ++ fileName
