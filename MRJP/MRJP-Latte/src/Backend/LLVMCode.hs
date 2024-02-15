module Backend.LLVMCode where

import qualified Data.Bifunctor
import Data.Char (ord)
import Data.List (intercalate)
import Debug.Trace (trace)
import qualified Frontend.AbstractSyntaxTree as AST
import Numeric (showHex)

spacing :: Int
spacing = 4

class Codegen a where
  codegen :: a -> String

class PhiReplacer a where
  replacePhi :: (String, String) -> a -> a

data Type
  = Void
  | Int
  | Bool
  | PtrAsInt
  | Char
  | String
  | Class String
  | Array Type
  | Ptr Type
  | Tuple [Type]
  | UnsafeRaw String
  | Fun Type [Type]
  deriving (Eq, Ord, Show)

voidPtrType :: Type
voidPtrType = Ptr Char

fromAstType :: AST.ASTType -> Type
fromAstType AST.ASTInt = Int
fromAstType AST.ASTBool = Bool
fromAstType AST.ASTStr = String
fromAstType (AST.ASTClass name) = Class name
fromAstType (AST.ASTArr astType) = Array $ fromAstType astType
fromAstType AST.ASTVoid = Void

mapIntercalated :: (a -> String) -> [a] -> String
mapIntercalated f = intercalate ", " . map f

mapIntercalatedCodegen :: Codegen a => [a] -> String
mapIntercalatedCodegen = mapIntercalated codegen

data Program = Program
  { progClasses :: [Class],
    progFunctions :: [Function],
    progStrings :: [StringLiteral]
  }
  deriving (Eq, Ord, Show)

data StringLiteral
  = StringLiteral String String
  deriving (Eq, Ord, Show)

data Class = ClassObj
  { className :: String,
    classFields :: [(Type, String)],
    classVTableType :: [Type],
    classVTable :: [String]
  }
  deriving (Eq, Ord, Show)

data Function = Function
  { funName :: String,
    funArgs :: [(String, Type)],
    funRet :: Type,
    funBlock :: [Instructions]
  }
  deriving (Eq, Ord, Show)

data Cond
  = Eq
  | Ne
  | Slt
  | Sle
  | Sgt
  | Sge
  deriving (Eq, Ord, Show)

data Instructions
  = Add String Type String String
  | Sub String Type String String
  | Mul String Type String String
  | Div String Type String String
  | Mod String Type String String
  | Neg String Type String
  | And String Type String String
  | Or String Type String String
  | Not String Type String
  | Icmp String Cond Type String String
  | Br String
  | Cbr String String String
  | Gep String Type String [String]
  | Phi String Type [(String, String)]
  | Label String
  | Call String Type String [(Type, String)]
  | CallM String Type String [(Type, String)]
  | Ret Type String
  | RetVoid
  | Comment String
  | Bitcast String String String String
  | BitcastT String Type String Type
  | Load String Type String
  | Store Type String Type String
  | Extractvalue String Type String [Int]
  | Insertvalue String Type String Type String [Int]
  | Ptrtoint String Type String Type
  deriving (Eq, Ord, Show)

instance Codegen Type where
  codegen Void = "void"
  codegen Int = "i32"
  codegen PtrAsInt = "i64"
  codegen Char = "i8"
  codegen Bool = "i1"
  codegen String = "%struct.String*"
  codegen (Class name) = "%" ++ name ++ "*"
  codegen (Array t) = "{ i32, " ++ codegen t ++ "* }*"
  codegen (Ptr t) = codegen t ++ "*"
  codegen (Tuple ts) = "{ " ++ mapIntercalatedCodegen ts ++ " }"
  codegen (UnsafeRaw s) = s
  codegen (Fun ret args) = codegen ret ++ "(" ++ mapIntercalatedCodegen args ++ ")*"

instance Codegen Program where
  codegen (Program classes functions strings) =
    unlines $
      map codegen strings ++ ["\n"]
        ++ map codegenLinted classes
        ++ ["\n"]
        ++ map codegenLinted functions

codegenLinted :: Codegen a => a -> String
codegenLinted = (++ "\n") . codegen

stringLiteralAllowedChars :: String
stringLiteralAllowedChars = ['a' .. 'z'] ++ ['A' .. 'Z'] ++ ['0' .. '9'] ++ [' ']

instance Codegen StringLiteral where
  codegen (StringLiteral name value) =
    name ++ " = internal constant [" ++ show (length value) ++ " x i8] c\"" ++ concatMap toHex value ++ "\""
    where
      toHex :: Char -> String
      toHex a | a `elem` stringLiteralAllowedChars = [a]
      toHex x = "\\" ++ getHex (ord x)
      getHex :: Int -> String
      getHex x = do
        let result = showHex x ""
        if length result == 1 then "0" ++ result else result

instance Codegen Class where
  codegen (ClassObj name fields vtableTypes vtable) =
    let vTableName = "%" ++ name ++ ".vtable" in
    let vTableData = "@" ++ name ++ ".vtable.data" in
    "%" ++ name ++ " = type { " ++ mapIntercalatedCodegen (Ptr (UnsafeRaw vTableName) : map fst fields) ++ " }\n"
      ++ vTableName ++ " = type { " ++ mapIntercalatedCodegen vtableTypes ++ " }\n\n"
      ++ vTableData ++ " = global " ++ vTableName ++ " {" ++ mapIntercalated makeVMet (zip vtableTypes vtable) ++ "\n}"
    where
      makeVMet :: (Type, String) -> String
      makeVMet (vTableName, vtable) = "\n    " ++ codegen vTableName ++ " @" ++ vtable

instance Codegen Function where
  codegen (Function name args ret block) =
    "define " ++ codegen ret ++ " @" ++ name ++ "(" ++ codegenArgs args ++ ") {\n"
      ++ unlines (map codegenAndIndent block)
      ++ "}"
    where
      codegenAndIndent s@Label {} = codegen s
      codegenAndIndent s@Comment {} = (replicate (spacing `div` 2) ' ' ++) $ codegen s
      codegenAndIndent s = (replicate spacing ' ' ++) $ codegen s

codegenArgs :: [(String, Type)] -> String
codegenArgs = mapIntercalated (\(name, t) -> codegen t ++ " %" ++ name)

instance Codegen Instructions where
  codegen (Add r t a b) = r ++ " = add " ++ codegen t ++ " " ++ a ++ ", " ++ b
  codegen (Sub r t a b) = r ++ " = sub " ++ codegen t ++ " " ++ a ++ ", " ++ b
  codegen (Mul r t a b) = r ++ " = mul " ++ codegen t ++ " " ++ a ++ ", " ++ b
  codegen (Div r t a b) = r ++ " = sdiv " ++ codegen t ++ " " ++ a ++ ", " ++ b
  codegen (Mod r t a b) = r ++ " = srem " ++ codegen t ++ " " ++ a ++ ", " ++ b
  codegen (Neg r t a) = r ++ " = sub " ++ codegen t ++ " 0, " ++ a
  codegen (And r t a b) = r ++ " = and " ++ codegen t ++ " " ++ a ++ ", " ++ b
  codegen (Or r t a b) = r ++ " = or " ++ codegen t ++ " " ++ a ++ ", " ++ b
  codegen (Not r t a) = r ++ " = " ++ "xor " ++ codegen t ++ " " ++ a ++ ", 1"
  codegen (Icmp r cond t a b) = r ++ " = icmp " ++ codegen cond ++ " " ++ codegen t ++ " " ++ a ++ ", " ++ b
  codegen (Br label) = "br label %" ++ label
  codegen (Cbr cond trueLabel falseLabel) = "br i1 " ++ cond ++ ", label %" ++ trueLabel ++ ", label %" ++ falseLabel
  codegen (Gep r t ptr idx) = r ++ " = getelementptr " ++ codegen t ++ ", " ++ codegen t ++ "* " ++ ptr ++ ", " ++ mapIntercalated (\x -> codegen Int ++ " " ++ x) idx
  codegen (Phi r t values) = r ++ " = phi " ++ codegen t ++ " " ++ intercalate ", " (map (\(label, val) -> "[" ++ val ++ ", %" ++ label ++ "]") values)
  codegen (Label label) = label ++ ":"
  codegen (Call _ Void name args) = "call " ++ codegen Void ++ " @" ++ name ++ "(" ++ intercalate ", " (map (\(ty, op) -> codegen ty ++ " " ++ op) args) ++ ")"
  codegen (CallM _ Void name args) = "call " ++ codegen Void ++ " " ++ name ++ "(" ++ intercalate ", " (map (\(ty, op) -> codegen ty ++ " " ++ op) args) ++ ")"
  codegen (Call ret t name args) = ret ++ " = call " ++ codegen t ++ " @" ++ name ++ "(" ++ intercalate ", " (map (\(ty, op) -> codegen ty ++ " " ++ op) args) ++ ")"
  codegen (CallM ret t name args) = ret ++ " = call " ++ codegen t ++ " " ++ name ++ "(" ++ intercalate ", " (map (\(ty, op) -> codegen ty ++ " " ++ op) args) ++ ")"
  codegen (Ret t val) = "ret " ++ codegen t ++ " " ++ val
  codegen RetVoid = "ret void"
  codegen (Comment s) = "; " ++ s
  codegen (Bitcast r t a b) = r ++ " = bitcast " ++ t ++ " " ++ a ++ " to " ++ b
  codegen (BitcastT r t a b) = r ++ " = bitcast " ++ codegen t ++ " " ++ a ++ " to " ++ codegen b
  codegen (Load r t ptr) = r ++ " = load " ++ codegen t ++ ", " ++ codegen t ++ "* " ++ ptr
  codegen (Store t val t' ptr) = "store " ++ codegen t ++ " " ++ val ++ ", " ++ codegen t' ++ " " ++ ptr
  codegen (Extractvalue r t ptr idxs) = r ++ " = extractvalue " ++ codegen t ++ " " ++ ptr ++ ", " ++ mapIntercalated show idxs
  codegen (Insertvalue r t ptr t' val idxs) = r ++ " = insertvalue " ++ codegen t ++ " " ++ ptr ++ ", " ++ codegen t' ++ " " ++ val ++ ", " ++ mapIntercalated show idxs
  codegen (Ptrtoint r t ptr t') = r ++ " = ptrtoint " ++ codegen t ++ " " ++ ptr ++ " to " ++ codegen t'

-- >>> codegen $ Insertvalue "a" (Ptr Int) "b" "c" "0"
-- "a = insertvalue i32* b, i32* c, 0"

instance Codegen Cond where
  codegen Eq = "eq"
  codegen Ne = "ne"
  codegen Slt = "slt"
  codegen Sle = "sle"
  codegen Sgt = "sgt"
  codegen Sge = "sge"

instance PhiReplacer Program where
  replacePhi (old, new) (Program classes functions strings) =
    Program
      (map (replacePhi (old, new)) classes)
      (map (replacePhi (old, new)) functions)
      strings

instance PhiReplacer Class where
  replacePhi (old, new) c = c

instance PhiReplacer Function where
  replacePhi (old, new) (Function name args ret block) =
    Function name args ret $ map (replacePhi (old, new)) block

instance PhiReplacer StringLiteral where
  replacePhi (old, new) s = s

replace :: Eq a => a -> a -> a -> a
replace old new x | x == old = new
replace old new x = x

instance PhiReplacer Instructions where
  replacePhi (old, new) (Add r t a b) = Add r t (replace old new a) (replace old new b)
  replacePhi (old, new) (Sub r t a b) = Sub r t (replace old new a) (replace old new b)
  replacePhi (old, new) (Mul r t a b) = Mul r t (replace old new a) (replace old new b)
  replacePhi (old, new) (Div r t a b) = Div r t (replace old new a) (replace old new b)
  replacePhi (old, new) (Mod r t a b) = Mod r t (replace old new a) (replace old new b)
  replacePhi (old, new) (Neg r t a) = Neg r t (replace old new a)
  replacePhi (old, new) (And r t a b) = And r t (replace old new a) (replace old new b)
  replacePhi (old, new) (Or r t a b) = Or r t (replace old new a) (replace old new b)
  replacePhi (old, new) (Not r t a) = Not r t (replace old new a)
  replacePhi (old, new) (Icmp r cond t a b) = Icmp r cond t (replace old new a) (replace old new b)
  replacePhi (old, new) (Br label) = Br (replace old new label)
  replacePhi (old, new) (Cbr cond trueLabel falseLabel) = Cbr (replace old new cond) (replace old new trueLabel) (replace old new falseLabel)
  replacePhi (old, new) (Gep r t ptr idx) = Gep r t (replace old new ptr) (map (replace old new) idx)
  replacePhi (old, new) (Phi r t values) = Phi r t $ map (Data.Bifunctor.bimap (replace old new) (replace old new)) values
  replacePhi (old, new) (Label label) = Label (replace old new label)
  replacePhi (old, new) (Call ret t name args) = Call ret t name $ map (Data.Bifunctor.second (replace old new)) args
  replacePhi (old, new) (CallM ret t name args) = CallM ret t name $ map (Data.Bifunctor.second (replace old new)) args
  replacePhi (old, new) (Ret t val) = Ret t (replace old new val)
  replacePhi (old, new) RetVoid = RetVoid
  replacePhi (old, new) (Comment s) = Comment s
  replacePhi (old, new) (Bitcast r t a b) = Bitcast r t (replace old new a) (replace old new b)
  replacePhi (old, new) (BitcastT r t a b) = BitcastT r t (replace old new a) b
  replacePhi (old, new) (Load r t ptr) = Load r t (replace old new ptr)
  replacePhi (old, new) (Store t val t' ptr) = Store t (replace old new val) t' (replace old new ptr)
  replacePhi (old, new) (Extractvalue r t ptr idxs) = Extractvalue r t (replace old new ptr) idxs
  replacePhi (old, new) (Insertvalue r t ptr t' val idxs) = Insertvalue r t (replace old new ptr) t' (replace old new val) idxs
  replacePhi (old, new) (Ptrtoint r t ptr t') = Ptrtoint r t (replace old new ptr) t'
