module Frontend.Commons where

import Latte.Abs
import System.IO

type Pos = BNFC'Position

getPosLine :: Pos -> Int
getPosLine (Just (line, _)) = line
getPosLine Nothing = -1

getPosCol :: Pos -> Int
getPosCol (Just (_, col)) = col
getPosCol Nothing = -1

posToString :: Pos -> String
posToString (Just (line, col)) = "line " ++ show line ++ ", column " ++ show col
posToString Nothing = "unknown position"

putStrLnStderr :: String -> IO ()
putStrLnStderr = hPutStrLn stderr
