module Frontend.Types where

import Control.Monad.Except
import Control.Monad.State
import qualified Data.Map as Map
import Data.Tuple.Extra ((&&&))
import qualified Frontend.AbstractSyntaxTree as AST
import Frontend.Commons
import Frontend.TypeErrors (ErrHolder (InvalidArrayType))
import Latte.Abs

type EM a = IO (Either ErrHolder a)

type TM a = StateT TypeCheckerEnv (ExceptT ErrHolder IO) a

data EnvType
  = EnvInt
  | EnvString
  | EnvBoolean
  | EnvVoid
  | EnvClass String
  | EnvArray EnvType
  deriving (Eq, Ord, Show)

toLangString :: EnvType -> String
toLangString EnvInt = "int"
toLangString EnvString = "string"
toLangString EnvBoolean = "boolean"
toLangString EnvVoid = "void"
toLangString (EnvClass name) = name
toLangString (EnvArray t) = toLangString t ++ "[]"

mapBNFCType :: Type -> TM EnvType
mapBNFCType Int {} = pure EnvInt
mapBNFCType Str {} = pure EnvString
mapBNFCType Bool {} = pure EnvBoolean
mapBNFCType Void {} = pure EnvVoid
mapBNFCType (Class _ (MIdent _ (Ident name))) = pure $ EnvClass name
mapBNFCType (Array pos t) = do
  t' <- mapBNFCType t
  t'str <- toErrStringRaw t
  unless (canBeArrays t') $ throwError $ InvalidArrayType pos t'str
  pure $ EnvArray t'
mapBNFCType Fun {} = error "Function cannot be part of user environment"

toErrStringRaw :: Type -> TM String
toErrStringRaw t = do
  t' <- mapBNFCType t
  pure $ toLangString t'

canBeArrays :: EnvType -> Bool
canBeArrays EnvInt = True
canBeArrays EnvString = True
canBeArrays EnvBoolean = True
canBeArrays EnvVoid = False
canBeArrays EnvClass {} = True
canBeArrays EnvArray {} = True -- TODO: restore before submission

isArray :: EnvType -> Bool
isArray EnvArray {} = True
isArray _ = False

isClass :: EnvType -> Bool
isClass EnvClass {} = True
isClass _ = False

-- name, return_type, [arg_type, arg_name]
data EnvFunction = EnvFunction
  { envFunName :: String,
    envFunRetType :: EnvType,
    envFunArgs :: [(EnvType, String)],
    envFunPos :: Pos
  }
  deriving (Eq, Ord, Show)

type GenericEnv a = Map.Map String a

data EnvClassField = EnvClassFieldObject
  { envClassFieldName :: String,
    envClassFieldType :: EnvType,
    envClassFieldPos :: Pos
  }
  deriving (Eq, Ord, Show)

data EnvClass = EnvClassObject
  { envClassName :: String,
    envClassFields :: GenericEnv EnvClassField,
    envClassMethods :: GenericEnv EnvFunction,
    envClassSuperClass :: Maybe String,
    envClassPos :: Pos
  }
  deriving (Eq, Ord, Show)

data ReturnType
  = Always
  | Sometimes
  | Never
  deriving (Eq, Ord, Show)

data ValueInfo = ValueInfo
  { valueInfoType :: EnvType,
    isValueAssignable :: Bool
  }

data TypeCheckerEnv = TypeCheckerEnv
  { envGetFunctions :: GenericEnv EnvFunction,
    envGetClasses :: GenericEnv EnvClass,
    envGetCurrentVariables :: GenericEnv EnvType,
    envGetCurrentClass :: Maybe EnvClass,
    envGetExpectedReturnType :: Maybe EnvType,
    evaluatorVarEnv :: GenericEnv AST.ASTType,
    evaluatorFunRetEnv :: GenericEnv AST.ASTType,
    evaluatorCurrentClass :: Maybe AST.ASTClass,
    evaluatorClasses :: GenericEnv AST.ASTClass,
    evaluatorExpectedReturn :: Maybe AST.ASTType
  }

mapEnvType :: EnvType -> AST.ASTType
mapEnvType EnvInt = AST.ASTInt
mapEnvType EnvString = AST.ASTStr
mapEnvType EnvBoolean = AST.ASTBool
mapEnvType EnvVoid = AST.ASTVoid
mapEnvType (EnvClass name) = AST.ASTClass name
mapEnvType (EnvArray t) = AST.ASTArr $ mapEnvType t

emptyTCEnv :: TypeCheckerEnv
emptyTCEnv =
  TypeCheckerEnv
    stdlib
    Map.empty
    Map.empty
    Nothing
    Nothing
    Map.empty
    ( Map.map (mapEnvType . envFunRetType) stdlib
    )
    Nothing
    Map.empty
    Nothing

stdlib :: GenericEnv EnvFunction
stdlib = Map.fromList $ map (envFunName &&& id) stdlibList
  where
    stdlibList =
      [ EnvFunction "printInt" EnvVoid [(EnvInt, "x")] Nothing,
        EnvFunction "printString" EnvVoid [(EnvString, "x")] Nothing,
        EnvFunction "error" EnvVoid [] Nothing,
        EnvFunction "readInt" EnvInt [] Nothing,
        EnvFunction "readString" EnvString [] Nothing
      ]

localState :: (TypeCheckerEnv -> TypeCheckerEnv) -> TM a -> TM a
localState envModifier localComputation = do
  savedState <- get
  put (envModifier savedState)
  localComputationResult <- localComputation
  put savedState
  pure localComputationResult
