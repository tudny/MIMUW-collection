module Utils where

import Data.Tuple (swap)
import System.Exit (ExitCode (..), exitWith)

exitCodeOfCode :: Int -> ExitCode
exitCodeOfCode 0 = ExitSuccess
exitCodeOfCode n = ExitFailure n

exitWithCode :: Int -> IO a
exitWithCode = exitWith . exitCodeOfCode

separator :: Char
separator = '/'

breakEnd :: (a -> Bool) -> [a] -> ([a], [a])
breakEnd f = swap . both reverse . break f . reverse
  where
    both f' (a, b) = (f' a, f' b)

dropLast :: [a] -> [a]
dropLast = reverse . drop 1 . reverse

-- Break a/b/c.ext into (a/b, c, ext)
decodeFilename :: String -> (String, String, String)
decodeFilename filename =
  let (path, filename') = fileIntoPathAndFilename filename
      (name, ext) = fileIntoNameAndExtension filename'
   in (path, name, ext)

-- Break a/b/c.ext into (a/b, c.ext)
fileIntoPathAndFilename :: String -> (String, String)
fileIntoPathAndFilename filename =
  let (path, filename') = breakEnd (== separator) filename
   in (dropLast path, filename')

fileIntoNameAndExtension :: String -> (String, String)
fileIntoNameAndExtension filename =
  let (name, ext) = breakEnd (== '.') filename
   in (dropLast name, ext)
