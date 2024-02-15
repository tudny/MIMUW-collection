module Frontend.AbstractSyntaxTree where

import qualified Data.Map as Map

newtype AST
  = Program [Def]
  deriving (Eq, Ord, Show)

type Field = (String, ASTType, Expr)

data Def
  = Function
      { defFunctions :: FunctionObj
      }
  | Class
      { defClassName :: String,
        defClassMethods :: [FunctionObj],
        defClassFields :: [Field]
      }
  deriving (Eq, Ord, Show)

data FunctionObj = FunctionObj
  { funObjName :: String,
    funObjArgs :: [Argument],
    funObjRet :: ASTType,
    funObjBlock :: Block
  }
  deriving (Eq, Ord, Show)

data Argument
  = Argument String ASTType
  deriving (Eq, Ord, Show)

newtype Block
  = Block [Stmt]
  deriving (Eq, Ord, Show)

data Stmt
  = SBlock Block
  | SVarDecl String ASTType Expr
  | SAssign Expr Expr
  | SExpr Expr
  | SIncr Expr
  | SDecr Expr
  | SReturn ASTType Expr
  | SReturnVoid
  | SIf Expr Block Block
  | SWhile Expr Block
  | SFor ASTType String Expr Block
  deriving (Eq, Ord, Show)

data Expr = Expr
  { exprType :: ASTType,
    exprValue :: Expr'
  }
  deriving (Eq, Ord, Show)

data Expr'
  = EVariable String
  | EInt Integer
  | EString String
  | EBool Bool
  | EApplication String [Expr]
  | EArrGet Expr Expr
  | EFieldGet Expr String
  | EArrLength ASTType Expr
  | EMethodCall Expr String [Expr]
  | ENewObject String
  | ENewArray Expr
  | ENull
  | ENeg Expr
  | ENot Expr
  | EMul Expr Expr
  | EDiv Expr Expr
  | EMod Expr Expr
  | EAdd Expr Expr
  | ESub Expr Expr
  | ELt Expr Expr
  | EGt Expr Expr
  | ELe Expr Expr
  | EGe Expr Expr
  | EEq Expr Expr
  | ENe Expr Expr
  | EAnd Expr Expr
  | EOr Expr Expr
  deriving (Eq, Ord, Show)

data ASTType
  = ASTInt
  | ASTStr
  | ASTBool
  | ASTVoid
  | ASTArr ASTType
  | ASTClass String
  deriving (Eq, Ord, Show)

data ASTFun = ASTFunObj
  { funName :: String,
    funArgs :: [ASTType],
    funRet :: ASTType
  }
  deriving (Eq, Ord, Show)

data ASTMethod = ASTMethodObj
  { methodFun :: ASTFun,
    methodIImplement :: Bool,
    methodImplementor :: String
  }
  deriving (Eq, Ord, Show)

alterFunctional :: (ASTFun -> ASTFun) -> (String -> String) -> ASTMethod -> ASTMethod
alterFunctional f g ASTMethodObj {methodFun = methodFun, methodIImplement = methodIImplement, methodImplementor = methodImplementor} =
  ASTMethodObj {methodFun = f methodFun, methodIImplement = methodIImplement, methodImplementor = g methodImplementor}

alterFunctionalT :: (ASTFun -> ASTFun) -> (String -> String) -> [ASTMethod] -> [ASTMethod]
alterFunctionalT f g = map (alterFunctional f g)

data ASTField = ASTFieldObj
  { classFieldName :: String,
    classFieldType :: ASTType
  }
  deriving (Eq, Ord, Show)

data ASTClass = ASTClassObj
  { className :: String,
    classFields :: [ASTField],
    classMethods :: [ASTMethod],
    classSuper :: Maybe String
  }
  deriving (Eq, Ord, Show)

data ASTEnv = ASTEnv
  { astEnvClasses :: Map.Map String ASTClass,
    astEnvFunctions :: Map.Map String ASTFun
  }
  deriving (Eq, Ord, Show)
