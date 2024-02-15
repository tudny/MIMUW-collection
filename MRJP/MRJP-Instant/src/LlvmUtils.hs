module LlvmUtils where

import Data.Text (pack, replace, unpack)

data Code
  = IADD String String String
  | ISUB String String String
  | IMUL String String String
  | IDIV String String String
  | IPRINT String
  | IALLOCA String
  | ISTORE String String
  | ILOAD String String

toString :: Code -> String
toString (IADD x y z) = x ++ " = add i32 " ++ y ++ ", " ++ z
toString (ISUB x y z) = x ++ " = sub i32 " ++ y ++ ", " ++ z
toString (IMUL x y z) = x ++ " = mul i32 " ++ y ++ ", " ++ z
toString (IDIV x y z) = x ++ " = sdiv i32 " ++ y ++ ", " ++ z
toString (IPRINT x) = "call void @printInt(i32 " ++ x ++ ")"
toString (IALLOCA x) = x ++ " = alloca i32"
toString (ISTORE x y) = "store i32 " ++ x ++ ", i32* " ++ y
toString (ILOAD x y) = x ++ " = load i32, i32* " ++ y

llvmTemplate :: String
llvmTemplate =
  "\
  \declare void @printInt(i32)\n\
  \define i32 @main() {\n\
  \%CODE_HERE%  ret i32 0\n\
  \}\n\
  \"

llvmCodePattern :: String
llvmCodePattern = "%CODE_HERE%"

parseLlvmTemplate :: [Code] -> IO String
parseLlvmTemplate code = do
  let codeString = unlines $ map (("  " ++) . toString) code
  let content = replace (pack llvmCodePattern) (pack codeString) (pack llvmTemplate)
  pure $ unpack content

saveLlvmFile :: String -> String -> String -> IO ()
saveLlvmFile path name content = do
  let filename = makeLlvmFileName path name
  writeFile filename content
  putStrLn $ "Generated " ++ filename

makeLlvmFileName :: String -> String -> String
makeLlvmFileName "" name = name ++ ".ll"
makeLlvmFileName path name = path ++ "/" ++ name ++ ".ll"
