module Frontend.TypeEvaluator where

import Control.Monad.Extra
import Control.Monad.State
import Data.List (sortBy)
import Data.List.Extra (find)
import qualified Data.Map as Map
import Data.Maybe (mapMaybe)
import Data.Tuple.Extra
import Frontend.AbstractSyntaxTree (AST (..), ASTClass (..), ASTEnv (..), ASTField (..), ASTFun (..))
import qualified Frontend.AbstractSyntaxTree as AST
import Frontend.ClassGraph
import Frontend.Commons (Pos)
import Frontend.Types (EnvClass (..), EnvClassField (EnvClassFieldObject), EnvFunction (..), EnvType (..), GenericEnv, TM, TypeCheckerEnv (..), localState, mapBNFCType, mapEnvType)
import qualified Latte.Abs as Abs

evalType :: Abs.Program -> TM (AST.AST, AST.ASTEnv)
evalType p = do
  env <- makeEnv
  let funs = Map.map funRet $ astEnvFunctions env
  let classes = astEnvClasses env
  ast <- withLocalFuns funs $ withClasses classes $ makeAst p
  pure (ast, env)

makeEnv :: TM AST.ASTEnv
makeEnv = do
  functions <- gets (mapFunctions . envGetFunctions)
  classes <- gets (mapClasses . envGetClasses)
  let classes' = enrichClassesWithSuperClasses classes
  pure $ AST.ASTEnv classes' functions

mapBnfcType :: Abs.Type -> TM AST.ASTType
mapBnfcType t = do
  envT <- mapBNFCType t
  pure $ mapEnvType envT

mapFunctions :: GenericEnv EnvFunction -> Map.Map String AST.ASTFun
mapFunctions =
  Map.map mapSingleFunction
  where
    mapSingleFunction :: EnvFunction -> AST.ASTFun
    mapSingleFunction (EnvFunction name retType args _) =
      AST.ASTFunObj name (map (mapEnvType . fst) args) (mapEnvType retType)

mapClasses :: GenericEnv EnvClass -> Map.Map String AST.ASTClass
mapClasses =
  Map.map mapSingleClass
  where
    mapSingleClass :: EnvClass -> AST.ASTClass
    mapSingleClass (EnvClassObject name fields methods superClassName _) =
      AST.ASTClassObj
        name
        (map mapSingleField $ Map.elems fields)
        (map (mapSingleMethod name) $ Map.elems methods)
        superClassName
    mapSingleField :: EnvClassField -> AST.ASTField
    mapSingleField (EnvClassFieldObject name t _) =
      AST.ASTFieldObj name (mapEnvType t)
    mapSingleMethod :: String -> EnvFunction -> AST.ASTMethod
    mapSingleMethod className (EnvFunction name retType args _) =
      AST.ASTMethodObj (AST.ASTFunObj name (map (mapEnvType . fst) args) (mapEnvType retType)) True className

-- monad utils --

addLocalVars :: [(String, AST.ASTType)] -> TM ()
addLocalVars vars =
  modify $ \env ->
    env
      { evaluatorVarEnv = Map.union (Map.fromList vars) (evaluatorVarEnv env)
      }

withLocalVars :: [(String, AST.ASTType)] -> TM a -> TM a
withLocalVars vars = do
  localState $ \env ->
    env
      { evaluatorVarEnv = Map.union (Map.fromList vars) (evaluatorVarEnv env)
      }

withLocalFuns :: Map.Map String AST.ASTType -> TM a -> TM a
withLocalFuns funs = do
  localState $ \env ->
    env
      { evaluatorFunRetEnv = Map.union funs (evaluatorFunRetEnv env)
      }

withClass :: AST.ASTClass -> TM a -> TM a
withClass c = do
  localState $ \env ->
    env
      { evaluatorCurrentClass = Just c
      }

withClasses :: Map.Map String AST.ASTClass -> TM a -> TM a
withClasses classes = do
  localState $ \env ->
    env
      { evaluatorClasses = Map.union classes (evaluatorClasses env)
      }

withExpectedReturnTy :: AST.ASTType -> TM a -> TM a
withExpectedReturnTy t = do
  localState $ \env ->
    env
      { evaluatorExpectedReturn = Just t
      }

-- defs --

defaultValueOfType :: AST.ASTType -> AST.Expr
defaultValueOfType AST.ASTInt = AST.Expr AST.ASTInt $ AST.EInt 0
defaultValueOfType AST.ASTStr = AST.Expr AST.ASTStr $ AST.EString ""
defaultValueOfType AST.ASTBool = AST.Expr AST.ASTBool $ AST.EBool False
defaultValueOfType ty@(AST.ASTArr t) = AST.Expr ty AST.ENull
defaultValueOfType ty@(AST.ASTClass name) = AST.Expr ty AST.ENull
defaultValueOfType AST.ASTVoid = error "void type has no default value"

-- ast makers --

tripleToPair :: (a, b, c) -> (a, b)
tripleToPair (a, b, _) = (a, b)

makeAst :: Abs.Program -> TM AST.AST
makeAst (Abs.Program _ topDefs) = do
  astTopDefs <- mapM makeTopDef topDefs
  pure $ AST.Program astTopDefs

makeTopDef :: Abs.TopDef -> TM AST.Def
makeTopDef (Abs.TopFnDef _ fnDef) = do
  fnObj <- makeFnDef fnDef
  pure $ AST.Function fnObj
makeTopDef (Abs.ClassExt _ (Abs.MIdent _ (Abs.Ident name)) _ items) = do
  Just currentClass@AST.ASTClassObj {classFields = classFields, classMethods = classMethods} <- gets (Map.lookup name . evaluatorClasses)
  let allFields = map (classFieldName &&& classFieldType) classFields
  let allMethods = Map.fromList $ map ((funName &&& funRet) . AST.methodFun) classMethods
  withClass currentClass $
    withLocalVars (("self", AST.ASTClass name) : allFields) $
      withLocalFuns allMethods $ do
        fields <- concatMapM makeClassField items
        methods <- concatMapM makeClassMethod items
        pure $ AST.Class name methods fields
makeTopDef (Abs.ClassDef p1 (Abs.MIdent p2 (Abs.Ident name)) items) =
  makeTopDef (Abs.ClassExt p1 (Abs.MIdent p2 (Abs.Ident name)) (Abs.MIdent p2 (Abs.Ident "")) items)

makeClassField :: Abs.ClassItem -> TM [(String, AST.ASTType, AST.Expr)]
makeClassField (Abs.ClassFieldNull _ t items) = do
  ty <- mapBnfcType t
  let defValue = defaultValueOfType ty
  mapM (makeClassFieldItem defValue ty) items
  where
    makeClassFieldItem :: AST.Expr -> AST.ASTType -> Abs.ClassField -> TM (String, AST.ASTType, AST.Expr)
    makeClassFieldItem defValue ty (Abs.ClassField _ (Abs.MIdent _ (Abs.Ident name))) = do
      pure (name, ty, defValue)
makeClassField Abs.ClassMethod {} = pure []

makeClassMethod :: Abs.ClassItem -> TM [AST.FunctionObj]
makeClassMethod (Abs.ClassMethod _ fnDef) = do
  fnDef' <- makeFnDef fnDef
  pure [fnDef']
makeClassMethod Abs.ClassFieldNull {} = pure []

makeFnDef :: Abs.FnDef -> TM AST.FunctionObj
makeFnDef (Abs.FnDef _ retType (Abs.MIdent _ (Abs.Ident name)) args block) = do
  astArgs <- mapM makeArg args
  let argsList = map (\(AST.Argument name t) -> (name, t)) astArgs
  retType' <- mapBNFCType retType
  astBlock <- withExpectedReturnTy (mapEnvType retType') $ withLocalVars argsList $ makeBlock block
  type' <- mapBnfcType retType
  pure $ AST.FunctionObj name astArgs type' astBlock

makeArg :: Abs.Arg -> TM AST.Argument
makeArg (Abs.Arg _ t (Abs.MIdent _ (Abs.Ident name))) = do
  t' <- mapBnfcType t
  pure $ AST.Argument name t'

makeBlock :: Abs.Block -> TM AST.Block
makeBlock (Abs.Block _ stmts) = do
  astStmts <- withLocalVars [] $ concatMapM makeStmt stmts
  pure $ AST.Block astStmts

makeStmt :: Abs.Stmt -> TM [AST.Stmt]
makeStmt Abs.Empty {} = pure []
makeStmt (Abs.BStmt _ block) = do
  astBlock <- makeBlock block
  pure [AST.SBlock astBlock]
makeStmt (Abs.Decl _ t items) = do
  t' <- mapBnfcType t
  concatMapM (makeItem t') items
makeStmt (Abs.Ass _ lVal rVal) = do
  lVal' <- makeExpr lVal
  rVal' <- makeExpr rVal
  pure [AST.SAssign lVal' rVal']
makeStmt (Abs.Incr _ lVal) = do
  lVal' <- makeExpr lVal
  pure [AST.SIncr lVal']
makeStmt (Abs.Decr _ lVal) = do
  lVal' <- makeExpr lVal
  pure [AST.SDecr lVal']
makeStmt (Abs.Ret _ e) = do
  e' <- makeExpr e
  Just retTy <- gets evaluatorExpectedReturn
  pure [AST.SReturn retTy e']
makeStmt Abs.VRet {} = do
  pure [AST.SReturnVoid]
makeStmt (Abs.Cond _ condE blockS) = do
  blockS' <- ensureBlock blockS
  condE' <- makeExpr condE
  pure [AST.SIf condE' blockS' $ AST.Block []]
makeStmt (Abs.CondElse _ condE blockSTrue blockSFalse) = do
  blockSTrue' <- ensureBlock blockSTrue
  blockSFalse' <- ensureBlock blockSFalse
  condE' <- makeExpr condE
  pure [AST.SIf condE' blockSTrue' blockSFalse']
makeStmt (Abs.While _ condE blockS) = do
  blockS' <- ensureBlock blockS
  condE' <- makeExpr condE
  pure [AST.SWhile condE' blockS']
makeStmt (Abs.For _ t (Abs.MIdent _ (Abs.Ident varName)) iterE blockS) = do
  varType <- mapBnfcType t
  blockS' <- withLocalVars [(varName, varType)] $ ensureBlock blockS
  iterE' <- makeExpr iterE
  pure [AST.SFor varType varName iterE' blockS']
makeStmt (Abs.SExp _ e) = do
  e' <- makeExpr e
  pure [AST.SExpr e']

makeItem :: AST.ASTType -> Abs.Item -> TM [AST.Stmt]
makeItem t (Abs.NoInit _ (Abs.MIdent _ (Abs.Ident name))) = do
  addLocalVars [(name, t)]
  pure [AST.SVarDecl name t (defaultValueOfType t)]
makeItem t (Abs.Init _ (Abs.MIdent _ (Abs.Ident name)) expr) = do
  astExpr <- makeExpr expr
  addLocalVars [(name, t)]
  pure [AST.SVarDecl name t astExpr]

ensureBlock :: Abs.Stmt -> TM AST.Block
ensureBlock (Abs.BStmt _ block) = makeBlock block
ensureBlock anyOtherStmt = makeBlock $ Abs.Block Nothing [anyOtherStmt]

makeExpr :: Abs.Expr -> TM AST.Expr
makeExpr (Abs.EVar _ (Abs.MIdent _ (Abs.Ident name))) = do
  varType <- takeVarType name
  pure $ AST.Expr varType $ AST.EVariable name
makeExpr (Abs.ELitInt _ int) = do
  pure $ AST.Expr AST.ASTInt $ AST.EInt int
makeExpr (Abs.ELitTrue _) = do
  pure $ AST.Expr AST.ASTBool $ AST.EBool True
makeExpr (Abs.ELitFalse _) = do
  pure $ AST.Expr AST.ASTBool $ AST.EBool False
makeExpr (Abs.EApp _ (Abs.FunctionCall _ (Abs.MIdent _ (Abs.Ident name)) exprs)) = do
  exprs' <- mapM makeExpr exprs
  retType <- takeFunRetType name
  pure $ AST.Expr retType $ AST.EApplication name exprs'
makeExpr (Abs.EString _ str) = do
  pure $ AST.Expr AST.ASTStr $ AST.EString str
makeExpr (Abs.ENullCast _ t) = do
  t' <- mapBnfcType t
  pure $ AST.Expr t' AST.ENull
makeExpr (Abs.EArrGet _ arrE idxE) = do
  arrE' <- makeExpr arrE
  idxE' <- makeExpr idxE
  let AST.ASTArr arrType = AST.exprType arrE'
  pure $ AST.Expr arrType $ AST.EArrGet arrE' idxE'
makeExpr (Abs.EClassGet _ objE (Abs.MIdent _ (Abs.Ident fieldName))) = do
  objE'@AST.Expr {AST.exprType = exprType} <- makeExpr objE
  case exprType of
    AST.ASTClass className -> do
      fieldType <- takeClassFieldType className fieldName
      pure $ AST.Expr fieldType $ AST.EFieldGet objE' fieldName
    AST.ASTArr arrT -> do
      pure $ AST.Expr AST.ASTInt $ AST.EArrLength arrT objE'
    _ -> error "type mismatch"
makeExpr (Abs.EClassMet _ objE (Abs.FunctionCall _ (Abs.MIdent _ (Abs.Ident methodName)) exprs)) = do
  objE' <- makeExpr objE
  exprs' <- mapM makeExpr exprs
  let AST.ASTClass className = AST.exprType objE'
  retType <- takeClassMethodType className methodName
  pure $ AST.Expr retType $ AST.EMethodCall objE' methodName exprs'
makeExpr (Abs.EArr _ t expr) = do
  expr' <- makeExpr expr
  t' <- mapBnfcType t
  pure $ AST.Expr (AST.ASTArr t') $ AST.ENewArray expr'
makeExpr (Abs.EClass _ (Abs.MIdent _ (Abs.Ident className))) = do
  pure $ AST.Expr (AST.ASTClass className) $ AST.ENewObject className
makeExpr (Abs.Neg _ expr) = do
  expr' <- makeExpr expr
  pure $ AST.Expr AST.ASTInt $ AST.ENeg expr'
makeExpr (Abs.Not _ expr) = do
  expr' <- makeExpr expr
  pure $ AST.Expr AST.ASTBool $ AST.ENot expr'
makeExpr (Abs.EMul _ expr1 op expr2) = do
  expr1' <- makeExpr expr1
  expr2' <- makeExpr expr2
  let op' = case op of
        Abs.Times _ -> AST.EMul
        Abs.Div _ -> AST.EDiv
        Abs.Mod _ -> AST.EMod
  pure $ AST.Expr AST.ASTInt $ op' expr1' expr2'
makeExpr (Abs.EAdd _ expr1 op expr2) = do
  expr1'@AST.Expr {AST.exprType = t1} <- makeExpr expr1
  expr2'@AST.Expr {AST.exprType = t2} <- makeExpr expr2
  unless (t1 == t2) $ error "type mismatch"
  let op' = case op of
        Abs.Plus _ -> AST.EAdd
        Abs.Minus _ -> AST.ESub
  pure $ AST.Expr t1 $ op' expr1' expr2'
makeExpr (Abs.ERel _ expr1 op expr2) = do
  expr1' <- makeExpr expr1
  expr2' <- makeExpr expr2
  let op' = case op of
        Abs.LTH _ -> AST.ELt
        Abs.LE _ -> AST.ELe
        Abs.GTH _ -> AST.EGt
        Abs.GE _ -> AST.EGe
        Abs.EQU _ -> AST.EEq
        Abs.NE _ -> AST.ENe
  pure $ AST.Expr AST.ASTBool $ op' expr1' expr2'
makeExpr (Abs.EAnd _ expr1 expr2) = do
  expr1' <- makeExpr expr1
  expr2' <- makeExpr expr2
  pure $ AST.Expr AST.ASTBool $ AST.EAnd expr1' expr2'
makeExpr (Abs.EOr _ expr1 expr2) = do
  expr1' <- makeExpr expr1
  expr2' <- makeExpr expr2
  pure $ AST.Expr AST.ASTBool $ AST.EOr expr1' expr2'

takeVarType :: String -> TM AST.ASTType
takeVarType name = do
  Just vars <- gets (Map.lookup name . evaluatorVarEnv)
  pure vars

takeFunRetType :: String -> TM AST.ASTType
takeFunRetType name = do
  Just funs <- gets (Map.lookup name . evaluatorFunRetEnv)
  pure funs

takeClassFieldType :: String -> String -> TM AST.ASTType
takeClassFieldType className fieldName = do
  Just classes <- gets (Map.lookup className . evaluatorClasses)
  let AST.ASTClassObj _ fields _ _ = classes
  let Just (AST.ASTFieldObj _ fieldType) = listLookup classFieldName fieldName fields
  pure fieldType

takeClassMethodType :: String -> String -> TM AST.ASTType
takeClassMethodType className methodName = do
  Just classes <- gets (Map.lookup className . evaluatorClasses)
  let AST.ASTClassObj _ _ methods _ = classes
  let Just (AST.ASTFunObj _ _ retType) = listLookup funName methodName (map AST.methodFun methods)
  pure retType

listLookup :: Eq k => (a -> k) -> k -> [a] -> Maybe a
listLookup f k = find (\x -> f x == k)
