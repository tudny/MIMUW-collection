module Backend.Renaming where

import Control.Monad.Reader
import Control.Monad.State
import qualified Data.Map as Map
import Frontend.AbstractSyntaxTree

type RM a = StateT RMEnv (ReaderT ASTEnv IO) a

data RMEnv = RMEnv
  { -- Map Variable Name to Next available number
    rmEnvVarCounter :: Map.Map String Int,
    -- Map Variable Name to Number of Currently Used Variable
    rmEnvCurrentVarRenames :: Map.Map String String,
    rmEnvCurrentClass :: Maybe String
  }
  deriving (Eq, Ord, Show)

emptyRMEnv :: RMEnv
emptyRMEnv = RMEnv Map.empty Map.empty Nothing

clearCollisions :: RM ()
clearCollisions = modify $ \env -> env {rmEnvVarCounter = Map.empty, rmEnvCurrentVarRenames = Map.empty}

localVariableNumbering :: RM a -> RM a
localVariableNumbering action = do
  keepRenames <- gets rmEnvCurrentVarRenames
  res <- action
  modify $ \env -> env {rmEnvCurrentVarRenames = keepRenames}
  pure res

modifyGlobalVarCounter :: (Map.Map String Int -> Map.Map String Int) -> RM ()
modifyGlobalVarCounter f = modify $ \env -> env {rmEnvVarCounter = f $ rmEnvVarCounter env}

type StrMap a = Map.Map String a

rename :: AST -> ASTEnv -> IO (AST, ASTEnv)
rename ast env = do
  env' <- renameEnv env
  ast' <- runReaderT (evalStateT (renameAST ast) emptyRMEnv) env'
  return (ast', env')

getCurrentCount :: String -> RM Int
getCurrentCount name = do
  levels <- gets rmEnvVarCounter
  case Map.lookup name levels of
    Just level -> pure level
    Nothing -> do
      let level = -1
      modifyGlobalVarCounter $ Map.insert name level
      pure level

getLocalRenaming :: String -> RM (Maybe String)
getLocalRenaming name = do
  renames <- gets rmEnvCurrentVarRenames
  pure $ Map.lookup name renames

getNameWithLevelOrMember :: String -> RM Expr'
getNameWithLevelOrMember name = do
  localRenaming <- getLocalRenaming name
  case localRenaming of
    Just localRenaming -> pure $ EVariable localRenaming
    Nothing -> do
      Just currentClass <- gets rmEnvCurrentClass
      if name == "self" then pure $ EVariable "self"
      else do
        let fieldName = enrichFieldNameFull currentClass name
        pure $ EFieldGet (Expr (ASTClass currentClass) (EVariable "self")) fieldName

getNextNameForVariable :: String -> RM String
getNextNameForVariable name = do
  level <- increaseLevel name
  let newName = enrichLocalName name level
  modify $ \env -> env {rmEnvCurrentVarRenames = Map.insert name newName $ rmEnvCurrentVarRenames env}
  pure newName

increaseLevel :: String -> RM Int
increaseLevel name = do
  level <- (1 +) <$> getCurrentCount name
  modifyGlobalVarCounter $ Map.insert name level
  pure level

localState :: (RMEnv -> RMEnv) -> RM a -> RM a
localState envModifier localComputation = do
  savedState <- get
  put (envModifier savedState)
  localComputationResult <- localComputation
  put savedState
  pure localComputationResult

withLocalClass :: String -> RM a -> RM a
withLocalClass className = localState $ \env -> env {rmEnvCurrentClass = Just className}

-- Naming convention:
-- Class names will be renamed to: class.xxx
-- Field will be renamed to: class.field.xxx
-- Methods will be renamed to: class.method.xxx
-- Global functions will be renamed to: function.xxx
-- Local variables will be renamed to: function.local.xxx.<level>

enrichClassName :: String -> String
enrichClassName = ("class." ++)

enrichFieldName :: String -> String -> String
enrichFieldName className = (("class." ++ className ++ ".field.") ++)

enrichFieldNameFull :: String -> String -> String
enrichFieldNameFull fullClassName = ((fullClassName ++ ".field.") ++)

enrichMethodName :: String -> String -> String
enrichMethodName fullClassName = ((fullClassName ++ ".method.") ++)

enrichFunctionName :: String -> String
enrichFunctionName = id

enrichLocalName :: String -> Int -> String
enrichLocalName name level = name ++ "." ++ show level

enrichAstType :: ASTType -> ASTType
enrichAstType (ASTArr subType) = ASTArr $ enrichAstType subType
enrichAstType (ASTClass name) = ASTClass $ enrichClassName name
enrichAstType t = t

renameEnv :: ASTEnv -> IO ASTEnv
renameEnv (ASTEnv classes functions) = do
  pure $ ASTEnv (renameClasses classes) (renameFunction functions)

renameClasses :: StrMap ASTClass -> StrMap ASTClass
renameClasses classes = do
  let classObjs = Map.elems classes
  Map.fromList $ map mapClassObj classObjs
  where
    mapClassObj :: ASTClass -> (String, ASTClass)
    mapClassObj (ASTClassObj name fields methods super) = do
      let name' = enrichClassName name
      let fields' = map (mapField name') fields
      let methods' = alterFunctionalT (mapMethod name') enrichClassName methods
      let super' = enrichClassName <$> super
      (name', ASTClassObj name' fields' methods' super')
    mapField :: String -> ASTField -> ASTField
    mapField className (ASTFieldObj name t) = do
      let name' = enrichFieldNameFull className name
      let t' = enrichAstType t
      ASTFieldObj name' t'
    mapMethod :: String -> ASTFun -> ASTFun
    mapMethod className (ASTFunObj name args ret) = do
      let name' = enrichMethodName className name
      let args' = map enrichAstType args
      let ret' = enrichAstType ret
      ASTFunObj name' args' ret'

renameFunction :: StrMap ASTFun -> StrMap ASTFun
renameFunction functions = do
  let funObjs = Map.elems functions
  Map.fromList $ map mapFunObj funObjs
  where
    mapFunObj :: ASTFun -> (String, ASTFun)
    mapFunObj (ASTFunObj name args ret) = do
      let name' = enrichFunctionName name
      let args' = map enrichAstType args
      let ret' = enrichAstType ret
      (name', ASTFunObj name' args' ret')

-- ast ---

renameAST :: AST -> RM AST
renameAST (Program defs) = do
  defs' <- mapM renameDef defs
  pure $ Program defs'

renameDef :: Def -> RM Def
renameDef (Class name methods fields) = do
  let name' = enrichClassName name
  withLocalClass name' $ do
    fields' <- mapM renameField fields
    methods' <- mapM (renameFunctionObj enrichMethodName) methods
    pure $ Class name' methods' fields'
renameDef (Function (FunctionObj name args retT block)) = do
  clearCollisions
  args' <- mapM renameArg args
  let retT' = enrichAstType retT
  block' <- renameBlock block
  pure $ Function $ FunctionObj (enrichFunctionName name) args' retT' block'

renameFunctionObj :: (String -> String -> String) -> FunctionObj -> RM FunctionObj
renameFunctionObj enr (FunctionObj name args ret body) = do
  clearCollisions
  Just currentClass <- gets rmEnvCurrentClass
  let name' = enr currentClass name
  args' <- mapM renameArg args
  let ret' = enrichAstType ret
  body' <- renameBlock body
  pure $ FunctionObj name' args' ret' body'

renameFunctionName :: String -> RM (String, Bool)
renameFunctionName functionName = do
  let functionRenamed = enrichFunctionName functionName
  currentClass <- gets rmEnvCurrentClass
  case currentClass of
    Just currentClass -> do
      Just currentClassObj <- asks $ Map.lookup currentClass . astEnvClasses
      let methods = classMethods currentClassObj
      let methodRenamed = enrichMethodName currentClass functionName
      pure $
        if methodRenamed `elem` map (funName . methodFun) methods
          then (methodRenamed, True)
          else (functionRenamed, False)
    Nothing -> pure (functionRenamed, False)

renameArg :: Argument -> RM Argument
renameArg (Argument name ty) = do
  name' <- getNextNameForVariable name
  let ty' = enrichAstType ty
  pure $ Argument name' ty'

renameField :: Field -> RM Field
renameField (name, ty, exp) = do
  exp' <- renameExpr exp
  Just currentClass <- gets rmEnvCurrentClass
  pure (enrichFieldNameFull currentClass name, enrichAstType ty, exp')

renameBlock :: Block -> RM Block
renameBlock (Block stmts) = do
  localVariableNumbering $ do
    stmts' <- mapM renameStmt stmts
    pure $ Block stmts'

renameStmt :: Stmt -> RM Stmt
renameStmt (SBlock block) = do
  block' <- renameBlock block
  pure $ SBlock block'
renameStmt (SVarDecl name ty value) = do
  value' <- renameExpr value
  name' <- getNextNameForVariable name
  let ty' = enrichAstType ty
  pure $ SVarDecl name' ty' value'
renameStmt (SAssign lhs rhs) = do
  rhs' <- renameExpr rhs
  lhs' <- renameExpr lhs
  pure $ SAssign lhs' rhs'
renameStmt (SExpr exp) = do
  exp' <- renameExpr exp
  pure $ SExpr exp'
renameStmt (SIncr exp) = do
  exp' <- renameExpr exp
  pure $ SIncr exp'
renameStmt (SDecr exp) = do
  exp' <- renameExpr exp
  pure $ SDecr exp'
renameStmt (SReturn ty exp) = do
  exp' <- renameExpr exp
  let ty' = enrichAstType ty
  pure $ SReturn ty' exp'
renameStmt SReturnVoid = pure SReturnVoid
renameStmt (SIf cond block1 block2) = do
  cond' <- renameExpr cond
  block1' <- renameBlock block1
  block2' <- renameBlock block2
  pure $ SIf cond' block1' block2'
renameStmt (SWhile cond block) = do
  cond' <- renameExpr cond
  block' <- renameBlock block
  pure $ SWhile cond' block'
renameStmt (SFor ty name exp block) = do
  exp' <- renameExpr exp
  localVariableNumbering $ do
    name' <- getNextNameForVariable name
    block' <- renameBlock block
    pure $ SFor (enrichAstType ty) name' exp' block'

renameExpr :: Expr -> RM Expr
renameExpr (Expr ty expr) = do
  expr' <- renameExpr' expr
  pure $ Expr (enrichAstType ty) expr'

renameExpr' :: Expr' -> RM Expr'
renameExpr' (EVariable name) = do
  getNameWithLevelOrMember name
renameExpr' (EInt i) = pure $ EInt i
renameExpr' (EString s) = pure $ EString s
renameExpr' (EBool b) = pure $ EBool b
renameExpr' (EApplication name args) = do
  (name', isMethod) <- renameFunctionName name
  args' <- mapM renameExpr args
  if isMethod
    then do
      Just currentClass <- gets rmEnvCurrentClass
      pure $ EMethodCall (Expr (ASTClass currentClass) (EVariable "self")) name' args'
    else
      pure $ EApplication name' args'
renameExpr' (EArrGet arr idx) = do
  arr' <- renameExpr arr
  idx' <- renameExpr idx
  pure $ EArrGet arr' idx'
renameExpr' (EFieldGet obj name) = do
  obj'@Expr {exprType = ASTClass className} <- renameExpr obj
  pure $ EFieldGet obj' (enrichFieldNameFull className name)
renameExpr' (EArrLength ty arr) = do
  arr' <- renameExpr arr
  pure $ EArrLength (enrichAstType ty) arr'
renameExpr' (EMethodCall obj name args) = do
  obj'@Expr {exprType = ASTClass className} <- renameExpr obj
  args' <- mapM renameExpr args
  pure $ EMethodCall obj' (enrichMethodName className name) args'
renameExpr' (ENewObject name) = pure $ ENewObject (enrichClassName name)
renameExpr' (ENewArray exp) = do
  exp' <- renameExpr exp
  pure $ ENewArray exp'
renameExpr' ENull = pure ENull
renameExpr' (ENeg exp) = do
  exp' <- renameExpr exp
  pure $ ENeg exp'
renameExpr' (ENot exp) = do
  exp' <- renameExpr exp
  pure $ ENot exp'
renameExpr' (EMul lhs rhs) = do
  lhs' <- renameExpr lhs
  rhs' <- renameExpr rhs
  pure $ EMul lhs' rhs'
renameExpr' (EDiv lhs rhs) = do
  lhs' <- renameExpr lhs
  rhs' <- renameExpr rhs
  pure $ EDiv lhs' rhs'
renameExpr' (EMod lhs rhs) = do
  lhs' <- renameExpr lhs
  rhs' <- renameExpr rhs
  pure $ EMod lhs' rhs'
renameExpr' (EAdd lhs rhs) = do
  lhs' <- renameExpr lhs
  rhs' <- renameExpr rhs
  pure $ EAdd lhs' rhs'
renameExpr' (ESub lhs rhs) = do
  lhs' <- renameExpr lhs
  rhs' <- renameExpr rhs
  pure $ ESub lhs' rhs'
renameExpr' (ELt lhs rhs) = do
  lhs' <- renameExpr lhs
  rhs' <- renameExpr rhs
  pure $ ELt lhs' rhs'
renameExpr' (EGt lhs rhs) = do
  lhs' <- renameExpr lhs
  rhs' <- renameExpr rhs
  pure $ EGt lhs' rhs'
renameExpr' (ELe lhs rhs) = do
  lhs' <- renameExpr lhs
  rhs' <- renameExpr rhs
  pure $ ELe lhs' rhs'
renameExpr' (EGe lhs rhs) = do
  lhs' <- renameExpr lhs
  rhs' <- renameExpr rhs
  pure $ EGe lhs' rhs'
renameExpr' (EEq lhs rhs) = do
  lhs' <- renameExpr lhs
  rhs' <- renameExpr rhs
  pure $ EEq lhs' rhs'
renameExpr' (ENe lhs rhs) = do
  lhs' <- renameExpr lhs
  rhs' <- renameExpr rhs
  pure $ ENe lhs' rhs'
renameExpr' (EAnd lhs rhs) = do
  lhs' <- renameExpr lhs
  rhs' <- renameExpr rhs
  pure $ EAnd lhs' rhs'
renameExpr' (EOr lhs rhs) = do
  lhs' <- renameExpr lhs
  rhs' <- renameExpr rhs
  pure $ EOr lhs' rhs'
