module JvmUtils where

import Data.Text (pack, replace, unpack)

data Code
  = IPUSH Integer
  | IADD
  | ISUB
  | IMUL
  | IDIV
  | ISTORE Integer
  | ILOAD Integer
  | IRETURN
  | IPRINT
  | SWAP

toString :: Code -> String
toString (IPUSH (-1)) = "iconst_m1"
toString (IPUSH x) | x >= 0 && x <= 5 = "iconst_" ++ show x
toString (IPUSH x) | x >= -128 && x <= 127 = "bipush " ++ show x
toString (IPUSH x) | x >= -32768 && x <= 32767 = "sipush " ++ show x
toString (IPUSH i) = "ldc " ++ show i
toString IADD = "iadd"
toString ISUB = "isub"
toString IMUL = "imul"
toString IDIV = "idiv"
toString (ISTORE i) | i >= 0 && i <= 3 = "istore_" ++ show i
toString (ISTORE i) = "istore " ++ show i
toString (ILOAD i) | i >= 0 && i <= 3 = "iload_" ++ show i
toString (ILOAD i) = "iload " ++ show i
toString IRETURN = "ireturn"
toString IPRINT = "invokestatic Runtime/printInt(I)V"
toString SWAP = "swap"

jvmTemplate :: String
jvmTemplate =
  "\
  \.class %CLASSNAME%\n\
  \.super java/lang/Object\n\
  \\n\
  \.method public <init>()V\n\
  \   aload_0\n\
  \   invokespecial java/lang/Object/<init>()V\n\
  \   return\n\
  \.end method\n\
  \\n\
  \.method public static main([Ljava/lang/String;)V\n\
  \.limit stack %STACK%\n\
  \.limit locals %LOCALS%\n\
  \%CODE_HERE%  return\n\
  \.end method"

jvmClassPattern :: String
jvmClassPattern = "%CLASSNAME%"

jvmCodePattern :: String
jvmCodePattern = "%CODE_HERE%"

jvmStackPattern :: String
jvmStackPattern = "%STACK%"

jvmLocalsPattern :: String
jvmLocalsPattern = "%LOCALS%"

pareseJvmTemplate :: String -> Integer -> Integer -> [Code] -> IO String
pareseJvmTemplate classname stackSize localsSize code = do
  let template = jvmTemplate
  let classPlaced'' = replace (pack jvmClassPattern) (pack classname) (pack template)
  let classPlaced' = replace (pack jvmStackPattern) (pack $ show stackSize) classPlaced''
  let classPlaced = replace (pack jvmLocalsPattern) (pack $ show (localsSize + 1)) classPlaced'
  let stringifyCode = map (("  " ++) . toString) code
  let joinCode = unlines stringifyCode
  let textifyCode = pack joinCode
  let codePlaced = replace (pack jvmCodePattern) textifyCode classPlaced
  pure $ unpack codePlaced

saveJasminFile :: String -> String -> String -> IO ()
saveJasminFile path name content = do
  let jasminFile = makeJasminFileName path name
  writeFile jasminFile content
  putStrLn $ "Generated: " ++ jasminFile

makeJasminFileName :: String -> String -> String
makeJasminFileName "" name = name ++ ".j"
makeJasminFileName path name = path ++ "/" ++ name ++ ".j"
