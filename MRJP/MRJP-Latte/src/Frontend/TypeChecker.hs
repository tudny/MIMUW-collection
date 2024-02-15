module Frontend.TypeChecker where

import Control.Arrow (ArrowChoice (right))
import Control.Monad.Except
import Control.Monad.HT (liftJoin2)
import Control.Monad.State
import qualified Data.Map as Map
import Data.Maybe
import qualified Data.Set as Set
import Data.Tuple.Extra
import qualified Frontend.AbstractSyntaxTree as AST
import Frontend.Commons
import Frontend.TypeErrors
import Frontend.TypeEvaluator
import Frontend.Types
import Latte.Abs

typeChecker :: Program -> EM (AST.AST, AST.ASTEnv)
typeChecker program = do
  runExceptT $ evalStateT (validateProgram program) emptyTCEnv

-- internal --------------------------------------------------------------------

withProtectedState :: TM a -> TM a
withProtectedState = localState id

withLocalVariables :: Map.Map String EnvType -> TM a -> TM a
withLocalVariables variables =
  localState $ \env ->
    env
      { envGetCurrentVariables = Map.union variables (envGetCurrentVariables env)
      }

withLocalVariable :: String -> EnvType -> TM a -> TM a
withLocalVariable name t =
  localState $ \env ->
    env
      { envGetCurrentVariables = Map.insert name t (envGetCurrentVariables env)
      }

withLocalClass :: EnvClass -> TM a -> TM a
withLocalClass classObject =
  localState $ \env ->
    env
      { envGetCurrentClass = Just classObject
      }

withExpectedReturnType :: EnvType -> TM a -> TM a
withExpectedReturnType t =
  localState $ \env ->
    env
      { envGetExpectedReturnType = Just t
      }

checkExpectedType :: Pos -> EnvType -> TM ()
checkExpectedType pos t = do
  expectedType <- gets envGetExpectedReturnType
  case expectedType of
    Nothing -> pure ()
    Just expected -> do
      ofType <- t `isOfType` expected
      unless ofType $ throwError $ InvalidReturnType pos (toLangString t) (toLangString expected)

withLocalVariablesAddClass :: EnvClass -> Map.Map String EnvType -> TM a -> TM a
withLocalVariablesAddClass classObject variables =
  localState $ \env ->
    env
      { envGetCurrentClass = Just classObject,
        envGetCurrentVariables = Map.union variables (envGetCurrentVariables env)
      }

addLocalVariable :: String -> EnvType -> TM ()
addLocalVariable name t =
  modify $ \env ->
    env
      { envGetCurrentVariables = Map.insert name t (envGetCurrentVariables env)
      }

addLocalVariables :: [(String, EnvType)] -> TM ()
addLocalVariables variables =
  modify $ \env ->
    env
      { envGetCurrentVariables = Map.union (Map.fromList variables) (envGetCurrentVariables env)
      }

makeResult :: a -> EM a
makeResult = pure . Right

makeError :: ErrHolder -> EM a
makeError = pure . Left

getEnvFunction :: String -> TM (Maybe EnvFunction)
getEnvFunction name = do
  classEnv <- gets envGetCurrentClass
  case classEnv of
    Nothing -> findInEnv name
    Just classObject -> findInClass classObject name
  where
    findInEnv :: String -> TM (Maybe EnvFunction)
    findInEnv name = do
      env <- gets envGetFunctions
      pure $ Map.lookup name env
    findInClass :: EnvClass -> String -> TM (Maybe EnvFunction)
    findInClass classObject name = do
      maybeMethod <- findMethodInClass name (Just classObject)
      case maybeMethod of
        Nothing -> findInEnv name
        Just _ -> pure maybeMethod

findMethodInClass :: String -> Maybe EnvClass -> TM (Maybe EnvFunction)
findMethodInClass _ Nothing = pure Nothing
findMethodInClass name (Just EnvClassObject {envClassMethods = methods, envClassSuperClass = superClass}) = do
  case Map.lookup name methods of
    Just fun -> pure $ Just fun
    Nothing -> do
      maybeSuperClass <- fromMaybe Nothing <$> mapM getEnvClass superClass
      findMethodInClass name maybeSuperClass

findFieldInClass :: String -> Maybe EnvClass -> TM (Maybe EnvType)
findFieldInClass _ Nothing = pure Nothing
findFieldInClass name (Just EnvClassObject {envClassFields = fields, envClassSuperClass = superClass}) = do
  case Map.lookup name fields of
    Just field -> pure $ Just $ envClassFieldType field
    Nothing -> do
      maybeSuperClass <- fromMaybe Nothing <$> mapM getEnvClass superClass
      findFieldInClass name maybeSuperClass

getEnvVariable :: String -> TM (Maybe EnvType)
getEnvVariable name = do
  env <- gets envGetCurrentVariables
  case Map.lookup name env of
    Just t -> pure $ Just t
    Nothing -> gets envGetCurrentClass >>= findFieldInClass' name
  where
    findFieldInClass' :: String -> Maybe EnvClass -> TM (Maybe EnvType)
    findFieldInClass' _ Nothing = pure Nothing
    findFieldInClass' name' (Just EnvClassObject {envClassFields = fields, envClassSuperClass = superClass}) = do
      case Map.lookup name' fields of
        Just field -> pure $ Just $ envClassFieldType field
        Nothing -> do
          maybeSuperClass <- fromMaybe Nothing <$> mapM getEnvClass superClass
          findFieldInClass' name' maybeSuperClass

addEnvFunctions :: [EnvFunction] -> TM ()
addEnvFunctions functions =
  modify $ \env -> env {envGetFunctions = Map.union (envGetFunctions env) $ Map.fromList $ map (envFunName &&& id) functions}

getEnvClass :: String -> TM (Maybe EnvClass)
getEnvClass name = do
  env <- gets envGetClasses
  pure $ Map.lookup name env

putEnvClasses :: [EnvClass] -> TM ()
putEnvClasses classes =
  modify $ \env -> env {envGetClasses = Map.fromList $ map (envClassName &&& id) classes}

-- validators ------------------------------------------------------------------

validateProgram :: Program -> TM (AST.AST, AST.ASTEnv)
validateProgram p@(Program _ topDefs) = do
  collectFunctions topDefs
  collectClasses topDefs
  validateFields
  checkProgram p
  evalType p

validateFields :: TM ()
validateFields = do
  classes <- gets envGetClasses
  mapM_ validateClassFields classes
  where
    validateClassFields :: EnvClass -> TM ()
    validateClassFields EnvClassObject {envClassName = _, envClassFields = fields, envClassPos = _} = do
      mapM_ validateClassField $ Map.elems fields
      where
        validateClassField :: EnvClassField -> TM ()
        validateClassField EnvClassFieldObject {envClassFieldName = fieldName, envClassFieldType = fieldType, envClassFieldPos = fieldPos} = do
          case fieldType of
            EnvInt -> pure ()
            EnvString -> pure ()
            EnvBoolean -> pure ()
            EnvVoid -> throwError $ InvalidVoidType fieldPos
            EnvClass name -> do
              maybeClass <- getEnvClass name
              case maybeClass of
                Nothing -> throwError $ UnknownClass name fieldPos
                Just _ -> pure ()
            EnvArray t -> validateClassField EnvClassFieldObject {envClassFieldName = fieldName, envClassFieldType = t, envClassFieldPos = fieldPos}
          case fieldName of
            "self" -> throwError $ InvalidFieldName fieldPos fieldName
            _ -> pure ()

collectFunctions :: [TopDef] -> TM ()
collectFunctions topDefs = do
  allFunctions <- catMaybes <$> mapM collectFunction topDefs
  builtins <- gets (Map.elems . envGetFunctions)
  catchError (checkDuplicateNames DuplicateFunctionName (allFunctions ++ builtins) envFunName envFunPos) handler
  addEnvFunctions allFunctions
  checkIfHasMain ()
  where
    collectFunction :: TopDef -> TM (Maybe EnvFunction)
    collectFunction ClassDef {} = pure Nothing
    collectFunction ClassExt {} = pure Nothing
    collectFunction (TopFnDef _ (FnDef _ retType (MIdent pos (Ident name)) args _)) = do
      validateFunArgListDuplicates args
      retType' <- mapBNFCType retType
      args' <- mapM mapArg args
      pure $
        Just $
          EnvFunction
            { envFunName = name,
              envFunRetType = retType',
              envFunArgs = args',
              envFunPos = pos
            }
    mapArg :: Arg -> TM (EnvType, String)
    mapArg (Arg _ t (MIdent _ (Ident name))) = do
      t' <- mapBNFCType t
      pure (t', name)
    handler :: ErrHolder -> TM ()
    handler (DuplicateFunctionName name positions)
      | Nothing `elem` positions = throwError $ DuplicateBuiltinFunctionName name (filter (/= Nothing) positions)
    handler err = throwError err

checkDuplicateNames :: (String -> [Pos] -> ErrHolder) -> [a] -> (a -> String) -> (a -> Pos) -> TM ()
checkDuplicateNames errorConstructor objs toName toPos = do
  let collected = groupBy toName objs
  let collectedWithPoses = Map.map (map toPos) collected
  let collectedWithMany = Map.filter (\x -> length x > 1) collectedWithPoses
  mapM_ (\(name, positions) -> throwError $ errorConstructor name positions) $ Map.toList collectedWithMany

findElem :: Eq a => a -> [(a, b)] -> Maybe b
findElem _ [] = Nothing
findElem x ((y, z) : ys) = if x == y then Just z else findElem x ys

validateFunArgListDuplicates :: [Arg] -> TM ()
validateFunArgListDuplicates args = do
  foldM_ validateFunArgDuplicates [] args
  where
    validateFunArgDuplicates :: [(String, Pos)] -> Arg -> TM [(String, Pos)]
    validateFunArgDuplicates names (Arg _ _ (MIdent pos (Ident name))) = do
      case findElem name names of
        Just firstPos -> throwError $ DuplicateArgName firstPos pos name
        Nothing -> pure $ (name, pos) : names

checkIfHasMain :: () -> TM ()
checkIfHasMain () = do
  maybeMain <- getEnvFunction "main"
  case maybeMain of
    Nothing -> throwError NoMain
    Just EnvFunction {envFunRetType = EnvInt, envFunArgs = [], envFunPos = _} -> pure ()
    Just EnvFunction {envFunRetType = EnvInt, envFunArgs = _, envFunPos = pos} -> throwError $ InvalidMainArgs pos
    Just EnvFunction {envFunRetType = _, envFunArgs = _, envFunPos = pos} -> throwError $ InvalidMainRetType pos

assoccBy :: Ord k => (a -> k) -> [a] -> Map.Map k a
assoccBy f = Map.fromList . map (f &&& id)

collectClasses :: [TopDef] -> TM ()
collectClasses topDefs = do
  collectedClasses <- catMaybes <$> mapM collectClass topDefs
  validateClassNames collectedClasses
  validateCycles collectedClasses
  putEnvClasses collectedClasses
  validateParentFields
  validateParentMethods
  checkClassDuplicates
  where
    collectClass :: TopDef -> TM (Maybe EnvClass)
    collectClass TopFnDef {} = pure Nothing
    collectClass (ClassDef _ (MIdent pos (Ident className)) classItems) = do
      fields <- assoccBy envClassFieldName <$> collectClassFields classItems
      methods <- assoccBy envFunName <$> collectClassMethods classItems
      pure $ Just $ EnvClassObject className fields methods Nothing pos
    collectClass (ClassExt _ (MIdent pos (Ident className)) (MIdent _ (Ident baseClassName)) classItems) = do
      fields <- assoccBy envClassFieldName <$> collectClassFields classItems
      methods <- assoccBy envFunName <$> collectClassMethods classItems
      pure $ Just $ EnvClassObject className fields methods (Just baseClassName) pos
    validateClassNames :: [EnvClass] -> TM ()
    validateClassNames classes = do
      mapM_ validateClassName classes
    validateClassName :: EnvClass -> TM ()
    validateClassName EnvClassObject {envClassName = className, envClassPos = pos} = do
      case className of
        "self" -> throwError $ InvalidClassName pos className
        _ -> pure ()

collectClassFields :: [ClassItem] -> TM [EnvClassField]
collectClassFields classItems = do
  collectedFields <- concat <$> mapM collectClassField classItems
  checkClassFieldsDuplicates collectedFields
  pure collectedFields
  where
    collectClassField :: ClassItem -> TM [EnvClassField]
    collectClassField (ClassFieldNull _ t items) = do
      mapM mapItem items
      where
        mapItem :: ClassField -> TM EnvClassField
        mapItem (ClassField _ (MIdent pos (Ident name))) = do
          t' <- mapBNFCType t
          pure $ EnvClassFieldObject name t' pos
    collectClassField (ClassMethod _ (FnDef _ _ (MIdent _ (Ident _)) _ _)) = pure []

checkClassFieldsDuplicates :: [EnvClassField] -> TM ()
checkClassFieldsDuplicates classFields = do
  foldM_ checkClassFieldDuplicates [] classFields
  where
    checkClassFieldDuplicates :: [(String, Pos)] -> EnvClassField -> TM [(String, Pos)]
    checkClassFieldDuplicates names (EnvClassFieldObject name _ pos) = do
      case findElem name names of
        Just firstPos -> throwError $ DuplicateClassFieldName name firstPos pos
        Nothing -> pure $ (name, pos) : names

collectClassMethods :: [ClassItem] -> TM [EnvFunction]
collectClassMethods classItems = do
  collectedMethods <- catMaybes <$> mapM collectClassMethod classItems
  checkClassMethodsDuplicates collectedMethods
  pure collectedMethods
  where
    collectClassMethod :: ClassItem -> TM (Maybe EnvFunction)
    collectClassMethod ClassFieldNull {} = pure Nothing
    collectClassMethod (ClassMethod _ (FnDef _ retType (MIdent pos (Ident name)) args _)) = do
      validateFunArgListDuplicates args
      retType' <- mapBNFCType retType
      args' <- mapM mapArg args
      pure $
        Just $
          EnvFunction
            { envFunName = name,
              envFunRetType = retType',
              envFunArgs = args',
              envFunPos = pos
            }
    mapArg :: Arg -> TM (EnvType, String)
    mapArg (Arg _ t (MIdent pos (Ident name))) = do
      t' <- mapBNFCType t
      pure (t', name)

checkClassMethodsDuplicates :: [EnvFunction] -> TM ()
checkClassMethodsDuplicates classMethods = do
  foldM_ checkClassMethodDuplicates [] classMethods
  where
    checkClassMethodDuplicates :: [(String, Pos)] -> EnvFunction -> TM [(String, Pos)]
    checkClassMethodDuplicates names (EnvFunction name _ _ pos) = do
      case findElem name names of
        Just firstPos -> throwError $ DuplicateClassMethodName name firstPos pos
        Nothing -> pure $ (name, pos) : names

data VisitedState = NotVisited | Visiting | Visited deriving (Eq, Ord, Show)

type Visited = Map.Map String VisitedState

validateCycles :: [EnvClass] -> TM ()
validateCycles classes = do
  let classNames = map envClassName classes
  let classPos = Map.fromList $ zip classNames $ map envClassPos classes
  let classMap = Map.fromList $ zip classNames classes
  let visited = Map.fromList $ zip classNames (repeat NotVisited)
  foldM_ (validateCycle classMap classPos Nothing) visited classNames
  where
    validateCycle :: Map.Map String EnvClass -> Map.Map String Pos -> Pos -> Visited -> String -> TM Visited
    validateCycle classMap classPos parentPos visited className = do
      let currentPos = getClassPos classPos className
      case Map.lookup className visited of
        Just Visited -> pure visited
        Just Visiting -> throwError $ CycleInheritance className currentPos
        Just NotVisited -> do
          let newVisited = Map.insert className Visiting visited
          case Map.lookup className classMap of
            Nothing -> throwError $ UnknownParentClass className currentPos
            Just EnvClassObject {envClassSuperClass = Nothing} -> do
              pure $ Map.insert className Visited visited
            Just EnvClassObject {envClassSuperClass = Just baseClassName} -> do
              _ <- validateCycle classMap classPos currentPos newVisited baseClassName
              pure $ Map.insert className Visited visited
        Nothing -> throwError $ UnknownParentClass className parentPos
      where
        getClassPos :: Map.Map String Pos -> String -> Pos
        getClassPos classPos' className' =
          join $ Map.lookup className' classPos'

type FieldMap = Map.Map String (Map.Map String Pos)

validateParentFields :: TM ()
validateParentFields = do
  classes <- gets $ Map.elems . envGetClasses
  mapM_ checkSingleClass classes
  where
    checkSingleClass :: EnvClass -> TM ()
    checkSingleClass
      EnvClassObject
        { envClassName = className,
          envClassFields = classFields,
          envClassSuperClass = superClass
        } = do
        superClassFieldMap <- getSuperClassFieldSet superClass
        let currentClassFieldMap = Map.fromList $ map (envClassFieldName &&& envClassFieldPos) (Map.elems classFields)
        let commonFields = Set.intersection (Map.keysSet superClassFieldMap) (Map.keysSet currentClassFieldMap)
        unless (null commonFields) $ throwError $ DuplicateSuperClassFieldName className (fromJust superClass) (getFieldsWithPos currentClassFieldMap superClassFieldMap commonFields)
    getSuperClassFieldSet :: Maybe String -> TM (Map.Map String Pos)
    getSuperClassFieldSet Nothing = pure Map.empty
    getSuperClassFieldSet (Just superClassName) = do
      Just superClass <- getEnvClass superClassName
      let superSuperClass = envClassSuperClass superClass
      superClassFieldMap <- getSuperClassFieldSet superSuperClass
      pure $ Map.map envClassFieldPos (envClassFields superClass) `Map.union` superClassFieldMap
    getFieldWithPos :: Map.Map String Pos -> String -> Pos
    getFieldWithPos fieldMap name = fromJust $ Map.lookup name fieldMap
    getFieldsWithPos :: Map.Map String Pos -> Map.Map String Pos -> Set.Set String -> [(String, Pos, Pos)]
    getFieldsWithPos currentClassFieldMap superClassFieldMap commonFields =
      map (\name -> (name, getFieldWithPos currentClassFieldMap name, getFieldWithPos superClassFieldMap name)) $ Set.toList commonFields

checkClassDuplicates :: TM ()
checkClassDuplicates = do
  classes <- gets $ Map.elems . envGetClasses
  checkDuplicateNames DuplicateClassName classes envClassName envClassPos

validateParentMethods :: TM ()
validateParentMethods = do
  classes <- gets $ Map.elems . envGetClasses
  mapM_ checkSingleClass classes
  where
    checkSingleClass :: EnvClass -> TM ()
    checkSingleClass
      EnvClassObject
        { envClassName = className,
          envClassMethods = classMethods,
          envClassSuperClass = superClass
        } = do
        superClassMethodMap <- getSuperClassMethodSet superClass
        let currentClassMethodMap = Map.fromList $ map (envFunName &&& (envFunPos &&& (envFunRetType &&& (map fst . envFunArgs)))) (Map.elems classMethods)
        let commonMethods = Set.intersection (Map.keysSet superClassMethodMap) (Map.keysSet currentClassMethodMap)
        checkTypeOfCommon currentClassMethodMap superClassMethodMap commonMethods
    getSuperClassMethodSet :: Maybe String -> TM (Map.Map String (Pos, (EnvType, [EnvType])))
    getSuperClassMethodSet Nothing = pure Map.empty
    getSuperClassMethodSet (Just superClassName) = do
      Just superClass <- getEnvClass superClassName
      let superSuperClass = envClassSuperClass superClass
      superClassMethodMap <- getSuperClassMethodSet superSuperClass
      pure $ superClassMethodMap `Map.union` Map.map (envFunPos &&& (envFunRetType &&& (map fst . envFunArgs))) (envClassMethods superClass)
    getMethodWithPos :: Map.Map String Pos -> String -> Pos
    getMethodWithPos methodMap name = fromJust $ Map.lookup name methodMap
    getMethodsWithPos :: Map.Map String Pos -> Map.Map String Pos -> Set.Set String -> [(String, Pos, Pos)]
    getMethodsWithPos currentClassMethodMap superClassMethodMap commonMethods =
      map (\name -> (name, getMethodWithPos currentClassMethodMap name, getMethodWithPos superClassMethodMap name)) $ Set.toList commonMethods
    checkTypeOfCommon :: Map.Map String (Pos, (EnvType, [EnvType])) -> Map.Map String (Pos, (EnvType, [EnvType])) -> Set.Set String -> TM ()
    checkTypeOfCommon current super common = do
      let currentCut = Map.assocs $ Map.restrictKeys current common
      let superCut = Map.assocs $ Map.restrictKeys super common
      let zipped = zip currentCut superCut
      mapM_ checkSingleCommon zipped
    checkSingleCommon :: ((String, (Pos, (EnvType, [EnvType]))), (String, (Pos, (EnvType, [EnvType])))) -> TM ()
    checkSingleCommon ((name, (pos, (retType, args))), (superName, (superPos, (superRetType, superArgs)))) = do
      when (name /= superName) $ error "Names should be the same"
      unless (retType == superRetType) $ throwError $ InvalidSuperClassMethodRetType name pos (toLangString retType) superPos (toLangString superRetType)
      unless (args == superArgs) $ throwError $ InvalidSuperClassMethodArgs name pos (map toLangString args) superPos (map toLangString superArgs)

-- checkers --------------------------------------------------------------------

checkProgram :: Program -> TM ()
checkProgram (Program _ topDefs) = do
  mapM_ checkTopDef topDefs

checkTopDef :: TopDef -> TM ()
checkTopDef (TopFnDef _ (FnDef _ retType (MIdent pos (Ident name)) args block)) = do
  validateFunArgListDuplicates args
  validateFunArgList args
  expectedReturnType <- mapBNFCType retType
  argsMap <- collectArgs args
  withLocalVariables argsMap $ do
    (retType', retPos) <- withExpectedReturnType expectedReturnType $ validateFunBlock block
    when (expectedReturnType /= EnvVoid) $ validateReturnType retPos retType'
  where
    collectArgs :: [Arg] -> TM (Map.Map String EnvType)
    collectArgs args' = do
      args'' <- mapM mapArg args'
      pure $ Map.fromList args''
    mapArg :: Arg -> TM (String, EnvType)
    mapArg (Arg pos t (MIdent _ (Ident name'))) = do
      t' <- mapBNFCType t
      pure (name', t')
    validateReturnType :: Pos -> ReturnType -> TM ()
    validateReturnType pos Never = throwError $ MissingReturnStatementInFunctionBlock pos name
    validateReturnType pos Sometimes = throwError $ MissingReturnStatementInFunctionBlockSometimes pos name
    validateReturnType _ Always = pure ()
checkTopDef (ClassDef _ (MIdent pos (Ident className)) classItems) = do
  maybeClass <- getEnvClass className
  when (isNothing maybeClass) $ throwError $ UnknownClass className pos
  let Just classObj = maybeClass
  withLocalClass classObj $ do
    mapM_ (checkClassItem classObj) classItems
  where
    checkClassItem :: EnvClass -> ClassItem -> TM ()
    checkClassItem _ ClassFieldNull {} = pure ()
    checkClassItem classObj (ClassMethod _ fnDef) = do
      withLocalClass classObj $
        withLocalVariable "self" (EnvClass className) $
          checkTopDef (TopFnDef pos fnDef)
checkTopDef (ClassExt pos className _ classItems) = checkTopDef (ClassDef pos className classItems)

validateFunArgList :: [Arg] -> TM ()
validateFunArgList args = do
  mapM_ validateFunSingleArg args
  where
    validateFunSingleArg :: Arg -> TM ()
    validateFunSingleArg (Arg _ t _) = do
      envContainsType t

envContainsType :: Type -> TM ()
envContainsType (Int _) = pure ()
envContainsType (Str _) = pure ()
envContainsType (Bool _) = pure ()
envContainsType (Void pos) = throwError $ InvalidVoidType pos
envContainsType (Class _ (MIdent pos (Ident name))) = do
  maybeClass <- getEnvClass name
  case maybeClass of
    Nothing -> throwError $ UnknownClass name pos
    Just _ -> pure ()
envContainsType (Array _ t) = do
  envContainsType t
envContainsType (Fun pos _ _) = throwError $ InvalidFunType pos

classToEnvType :: Pos -> String -> TM ValueInfo
classToEnvType pos name = do
  maybeClass <- getEnvClass name
  case maybeClass of
    Nothing -> throwError $ UnknownClass name pos
    Just _ -> pureInfo (EnvClass name) False

pureInfo :: EnvType -> Bool -> TM ValueInfo
pureInfo t assignable = pure $ ValueInfo t assignable

firstNotNone :: Maybe a -> Maybe a -> Maybe a
firstNotNone Nothing Nothing = Nothing
firstNotNone (Just x) _ = Just x
firstNotNone _ (Just x) = Just x

joinTwoReturnTypesBranches :: (ReturnType, Pos) -> (ReturnType, Pos) -> TM (ReturnType, Pos)
joinTwoReturnTypesBranches (Never, pos1) (Never, pos2) = pure (Never, firstNotNone pos1 pos2)
joinTwoReturnTypesBranches (Never, pos1) (Sometimes, pos2) = pure (Sometimes, firstNotNone pos1 pos2)
joinTwoReturnTypesBranches (Never, pos1) (Always, pos2) = pure (Sometimes, firstNotNone pos1 pos2)
joinTwoReturnTypesBranches (Sometimes, pos1) (Never, pos2) = pure (Sometimes, firstNotNone pos1 pos2)
joinTwoReturnTypesBranches (Sometimes, pos1) (Sometimes, pos2) = pure (Sometimes, firstNotNone pos1 pos2)
joinTwoReturnTypesBranches (Sometimes, pos1) (Always, pos2) = pure (Sometimes, firstNotNone pos1 pos2)
joinTwoReturnTypesBranches (Always, pos1) (Never, pos2) = pure (Sometimes, firstNotNone pos2 pos1)
joinTwoReturnTypesBranches (Always, pos1) (Sometimes, pos2) = pure (Sometimes, firstNotNone pos2 pos1)
joinTwoReturnTypesBranches (Always, pos1) (Always, pos2) = pure (Always, firstNotNone pos1 pos2)

joinTwoReturnTypesSeq :: (ReturnType, Pos) -> (ReturnType, Pos) -> TM (ReturnType, Pos)
joinTwoReturnTypesSeq (Never, pos1) (Never, pos2) = pure (Never, firstNotNone pos1 pos2)
joinTwoReturnTypesSeq (Never, pos1) (Sometimes, pos2) = pure (Sometimes, firstNotNone pos1 pos2)
joinTwoReturnTypesSeq (Never, pos1) (Always, pos2) = pure (Always, firstNotNone pos1 pos2)
joinTwoReturnTypesSeq (Sometimes, pos1) (Never, pos2) = pure (Sometimes, firstNotNone pos1 pos2)
joinTwoReturnTypesSeq (Sometimes, pos1) (Sometimes, pos2) = pure (Sometimes, firstNotNone pos1 pos2)
joinTwoReturnTypesSeq (Sometimes, pos1) (Always, pos2) = pure (Always, firstNotNone pos1 pos2)
joinTwoReturnTypesSeq (Always, pos1) (Never, pos2) = pure (Always, firstNotNone pos2 pos1)
joinTwoReturnTypesSeq (Always, pos1) (Sometimes, pos2) = pure (Always, firstNotNone pos2 pos1)
joinTwoReturnTypesSeq (Always, pos1) (Always, pos2) = pure (Always, firstNotNone pos2 pos1)

validateFunBlock :: Block -> TM (ReturnType, Pos)
validateFunBlock (Block blockPos stmts) = do
  checkBlockLocalVars stmts
  checkTypesStmts blockPos stmts

checkBlockLocalVars :: [Stmt] -> TM ()
checkBlockLocalVars stmts = do
  collectedDecls <- concat <$> mapM collectLocalDecls stmts
  let mergedDecls = Map.map (map snd) $ groupBy fst collectedDecls
  void $ Map.traverseWithKey checkDeclDup mergedDecls
  checkSubBlocks stmts
  where
    collectLocalDecls :: Stmt -> TM [(String, Pos)]
    collectLocalDecls (Decl _ _ items) = mapM collectLocalDecl items
    collectLocalDecls _ = pure []
    collectLocalDecl :: Item -> TM (String, Pos)
    collectLocalDecl (NoInit _ (MIdent pos (Ident name))) = pure (name, pos)
    collectLocalDecl (Init _ (MIdent pos (Ident name)) _) = pure (name, pos)
    checkDeclDup :: String -> [Pos] -> TM ()
    checkDeclDup name ps@(_ : _ : _) = throwError $ DuplicateLocalVarName name ps
    checkDeclDup _ _ = pure ()
    checkSubBlocks :: [Stmt] -> TM ()
    checkSubBlocks = mapM_ checkSubBlock
    checkSubBlock :: Stmt -> TM ()
    checkSubBlock (BStmt _ (Block _ stmts')) = checkBlockLocalVars stmts'
    checkSubBlock _ = pure ()

groupBy :: Ord k => (v -> k) -> [v] -> Map.Map k [v]
groupBy key = Map.fromListWith (++) . map (key &&& pure)

checkTypesStmts :: Pos -> [Stmt] -> TM (ReturnType, Pos)
checkTypesStmts topPos stmts = do
  returnedTypes <- withProtectedState $ mapM checkTypes stmts
  findMismatches returnedTypes
  where
    findMismatches :: [(ReturnType, Pos)] -> TM (ReturnType, Pos)
    findMismatches [] = pure (Never, topPos)
    findMismatches [(t, pos)] = pure (t, pos)
    findMismatches ((t1, pos1) : (t2, pos2) : ts) = do
      joined <- joinTwoReturnTypesSeq (t1, pos1) (t2, pos2)
      findMismatches $ joined : ts

checkTypes :: Stmt -> TM (ReturnType, Pos)
checkTypes (Empty pos) = pure (Never, pos)
checkTypes (Ret pos e) = do
  _ <- calculateLiterals e
  ValueInfo retType _ <- checkTypeE e
  if retType == EnvVoid
    then throwError $ ReturnVoid pos
    else do
      checkExpectedType pos retType
      pure (Always, pos)
checkTypes (VRet pos) = do
  checkExpectedType pos EnvVoid
  pure (Always, pos)
checkTypes (BStmt _ (Block blockPos stmts)) = checkTypesStmts blockPos stmts
checkTypes (Decl pos t items) = do
  envContainsType t
  t' <- mapBNFCType t
  collectItemsToEnv t' items
  pure (Never, pos)
  where
    collectItemsToEnv :: EnvType -> [Item] -> TM ()
    collectItemsToEnv varType its = do
      mapM_ collectItem its
      where
        collectItem :: Item -> TM ()
        collectItem (NoInit _ (MIdent _ (Ident name))) = do
          addLocalVariable name varType
        collectItem (Init _ (MIdent _ (Ident name)) e) = do
          _ <- calculateLiterals e
          ValueInfo eType _ <- checkTypeE e
          ofType <- eType `isOfType` varType
          unless ofType $ throwError $ DeclTypeMismatch pos (toLangString eType) (toLangString varType)
          addLocalVariable name varType
checkTypes (Ass pos lVal expr) = do
  ValueInfo l assignable <- checkTypeE lVal
  _ <- calculateLiterals expr
  ValueInfo t _ <- checkTypeE expr
  unless assignable $ throwError $ CannotChangeNotLValue pos
  canAssignType <- t `isOfType` l
  unless canAssignType $ throwError $ CannotAssignType pos (toLangString l) (toLangString t)
  pure (Never, pos)
checkTypes (Incr pos lVal) = do
  ValueInfo l assignable <- checkTypeE lVal
  unless assignable $ throwError $ CannotChangeNotLValue pos
  if l /= EnvInt
    then throwError $ CannotIncrementNotInt pos (toLangString l)
    else pure (Never, pos)
checkTypes (Decr pos lVal) = do
  ValueInfo l assignable <- checkTypeE lVal
  unless assignable $ throwError $ CannotChangeNotLValue pos
  if l /= EnvInt
    then throwError $ CannotDecrementNotInt pos (toLangString l)
    else pure (Never, pos)
checkTypes (Cond pos cond body) = do
  preCondType <- calculateLiterals cond
  ValueInfo condType _ <- checkTypeE cond
  when (condType /= EnvBoolean) $ throwError $ IfCondShouldBeBoolean pos (toLangString condType)
  blockType <- withProtectedState $ checkTypes body
  pure $
    if preCondType == EvaluatedBool False
      then (Never, pos)
      else
        if preCondType == EvaluatedBool True
          then blockType
          else go blockType
  where
    go (Never, pos') = (Never, pos')
    go (Sometimes, pos') = (Sometimes, pos')
    go (Always, _) = (Sometimes, pos)
checkTypes (CondElse pos cond body elseBody) = do
  ValueInfo condType _ <- checkTypeE cond
  when (condType /= EnvBoolean) $ throwError $ IfCondShouldBeBoolean pos (toLangString condType)
  condVal <- calculateLiterals cond
  brt <- withProtectedState $ checkTypes body
  ert <- withProtectedState $ checkTypes elseBody
  go condVal brt ert
  where
    go (EvaluatedBool True) r _ = pure r
    go (EvaluatedBool False) _ r = pure r
    go _ (t1, p1) (t2, p2) = joinTwoReturnTypesBranches (t1, p1) (t2, p2)
checkTypes (While pos cond body) = do
  ValueInfo condType _ <- checkTypeE cond
  when (condType /= EnvBoolean) $ throwError $ WhileCondShouldBeBoolean pos (toLangString condType)
  condVal <- calculateLiterals cond
  brt <- withProtectedState $ checkTypes body
  pure $ go condVal brt
  where
    go (EvaluatedBool True) r = r
    go (EvaluatedBool False) _ = (Never, pos)
    go _ (Never, pos') = (Never, pos')
    go _ (Always, pos') = (Sometimes, pos')
    go _ (Sometimes, pos') = (Sometimes, pos')
checkTypes (SExp pos e) = calculateLiterals e >> checkTypeE e >> pure (Never, pos)
checkTypes (For pos varType (MIdent _ (Ident name)) iterableExpr body) = do
  ValueInfo iterableType _ <- checkTypeE iterableExpr
  unless (isArray iterableType) $ throwError $ ForIterableNotArray pos (toLangString iterableType)
  typeMapped <- mapBNFCType varType
  let EnvArray arrayType = iterableType
  ofType <- arrayType `isOfType` typeMapped
  unless ofType $ throwError $ ForIterableTypeMismatch pos (toLangString iterableType) (toLangString typeMapped)
  bodyType <- withLocalVariable name typeMapped $ checkTypes body
  go bodyType
  where
    go (Never, pos') = pure (Never, pos')
    go (Always, pos') = pure (Sometimes, pos')
    go (Sometimes, pos') = pure (Sometimes, pos')

data ExprValueOnLits
  = EvaluatedBool Bool
  | EvaluatedInt Integer
  | TooComplex
  deriving (Eq)

minIntBound :: Integer
minIntBound = -2147483648

maxIntBound :: Integer
maxIntBound = 2147483647

checkOverflow :: Pos -> Integer -> TM ()
checkOverflow pos value = do
  when (value < minIntBound || value > maxIntBound) $ throwError $ IntOverflow pos value

checkOverflowOnLits :: Pos -> ExprValueOnLits -> TM ExprValueOnLits
checkOverflowOnLits pos x@(EvaluatedInt i) = do
  checkOverflow pos i
  pure x
checkOverflowOnLits _ x = pure x

if' :: ExprValueOnLits -> TM ExprValueOnLits -> TM ExprValueOnLits -> TM ExprValueOnLits
if' (EvaluatedBool True) t _ = t
if' (EvaluatedBool False) _ f = f
if' _ _ _ = pure TooComplex

calculateLiterals :: Expr -> TM ExprValueOnLits
calculateLiterals (ELitInt pos i) = checkOverflowOnLits pos $ EvaluatedInt i
calculateLiterals (ELitTrue _) = pure $ EvaluatedBool True
calculateLiterals (ELitFalse _) = pure $ EvaluatedBool False
calculateLiterals (Neg pos e) = do
  v' <- go <$> calculateLiterals e
  checkOverflowOnLits pos v'
  where
    go (EvaluatedInt i) = EvaluatedInt $ negate i
    go _ = TooComplex
calculateLiterals (Not _ e) = go <$> calculateLiterals e
  where
    go (EvaluatedBool b) = EvaluatedBool $ not b
    go _ = TooComplex
calculateLiterals (EMul pos e1 op e2) = calculateTwoArgOpInt pos (checkZero op) (resolveMulOp op) e1 e2
calculateLiterals (EAdd pos e1 op e2) = calculateTwoArgOpInt pos False (resolveAddOp op) e1 e2
calculateLiterals (ERel _ e1 rel e2) = calculateTwoArgOpBool (resolveRelOp rel) e1 e2
calculateLiterals (EAnd _ e1 e2) = do
  val1 <- calculateLiterals e1
  if' val1 (calculateLiterals e2) (pure val1)
calculateLiterals (EOr _ e1 e2) = do
  val1 <- calculateLiterals e1
  if' val1 (pure val1) (calculateLiterals e2)
calculateLiterals _ = pure TooComplex

resolveMulOp :: MulOp -> (Integer -> Integer -> Integer)
resolveMulOp Times {} = (*)
resolveMulOp Div {} = div
resolveMulOp Mod {} = rem

checkZero :: MulOp -> Bool
checkZero Times {} = False
checkZero Div {} = True
checkZero Mod {} = True

resolveAddOp :: AddOp -> (Integer -> Integer -> Integer)
resolveAddOp Plus {} = (+)
resolveAddOp Minus {} = (-)

resolveRelOp :: RelOp -> (Integer -> Integer -> Bool)
resolveRelOp LTH {} = (<)
resolveRelOp LE {} = (<=)
resolveRelOp GTH {} = (>)
resolveRelOp GE {} = (>=)
resolveRelOp EQU {} = (==)
resolveRelOp NE {} = (/=)

calculateTwoArgOpInt :: Pos -> Bool -> (Integer -> Integer -> Integer) -> Expr -> Expr -> TM ExprValueOnLits
calculateTwoArgOpInt pos cz op e1 e2 = do
  v <- liftJoin2 (go cz) (calculateLiterals e1) (calculateLiterals e2)
  checkOverflowOnLits pos v
  where
    go :: Bool -> ExprValueOnLits -> ExprValueOnLits -> TM ExprValueOnLits
    go True _ (EvaluatedInt 0) = throwError $ DivisionByZero pos
    go _ (EvaluatedInt i1) (EvaluatedInt i2) = pure $ EvaluatedInt $ op i1 i2
    go _ _ _ = pure TooComplex

calculateTwoArgOpBool :: (Integer -> Integer -> Bool) -> Expr -> Expr -> TM ExprValueOnLits
calculateTwoArgOpBool op e1 e2 = go <$> calculateLiterals e1 <*> calculateLiterals e2
  where
    go (EvaluatedInt i1) (EvaluatedInt i2) = EvaluatedBool $ op i1 i2
    go _ _ = TooComplex

getVarType :: String -> TM (Maybe EnvType)
getVarType varName = do
  maybeVar <- getEnvVariable varName
  case maybeVar of
    Nothing -> pure Nothing
    Just t -> pure $ Just t

checkTypeE :: Expr -> TM ValueInfo
checkTypeE (EVar _ (MIdent namePos (Ident name))) = do
  currentClass <- gets envGetCurrentClass
  maybeType <- getVarType name
  case maybeType of
    Nothing -> throwError $ UnknownVar name namePos
    Just t -> pureInfo t (isNothing currentClass || name /= "self")
checkTypeE (ELitInt _ _) = pureInfo EnvInt False
checkTypeE (ELitTrue _) = pureInfo EnvBoolean False
checkTypeE (ELitFalse _) = pureInfo EnvBoolean False
checkTypeE (EApp _ (FunctionCall funPos (MIdent _ (Ident name)) args)) =
  validateFunctionCall name funPos args
checkTypeE (ENullCast pos t) = do
  envContainsType t
  castedType <- mapBNFCType t
  if isClass castedType
    then pureInfo (EnvClass $ toLangString castedType) False
    else
      if isArray castedType
        then pureInfo castedType False
        else throwError $ InvalidNullCast pos (toLangString castedType)
checkTypeE (EString _ _) = pureInfo EnvString False
checkTypeE (EArr pos t e) = do
  envContainsType t
  realT <- valueInfoType <$> checkTypeE e
  arrayType <- mapBNFCType t
  if realT == EnvInt
    then pureInfo (EnvArray arrayType) False
    else throwError $ InvalidArraySize pos (toLangString realT)
checkTypeE (EArrGet pos array index) = checkArrAccess pos array index
checkTypeE (EClass _ (MIdent pos (Ident name))) = classToEnvType pos name
checkTypeE (EClassGet pos lVal (MIdent _ (Ident name))) = getClassField pos lVal name
checkTypeE (EClassMet pos lVal (FunctionCall _ (MIdent _ (Ident funName)) args)) = do
  checkMethodCall pos lVal funName args
checkTypeE (Neg pos e) = do
  ValueInfo t _ <- checkTypeE e
  if t == EnvInt
    then pureInfo EnvInt False
    else throwError $ InvalidNegType pos (toLangString t)
checkTypeE (Not pos e) = do
  ValueInfo t _ <- checkTypeE e
  if t == EnvBoolean
    then pureInfo EnvBoolean False
    else throwError $ InvalidNotType pos (toLangString t)
checkTypeE (EMul _ e1 op e2) = do
  ValueInfo t1 _ <- checkTypeE e1
  ValueInfo t2 _ <- checkTypeE e2
  if t1 == EnvInt && t2 == EnvInt
    then pureInfo EnvInt False
    else throwError $ InvalidMulType (hasPosition op) (toLangString t1) (toLangString t2)
checkTypeE (EAdd _ e1 op e2) = do
  ValueInfo t1 _ <- checkTypeE e1
  ValueInfo t2 _ <- checkTypeE e2
  go op t1 t2 >>= (`pureInfo` False)
  where
    go :: AddOp -> EnvType -> EnvType -> TM EnvType
    go _ EnvInt EnvInt = pure EnvInt
    go Plus {} EnvString EnvString = pure EnvString
    go _ t1 t2 = throwError $ InvalidAddType (hasPosition op) (toLangString t1) (toLangString t2)
checkTypeE (ERel _ e1 op e2) = do
  ValueInfo t1 _ <- checkTypeE e1
  ValueInfo t2 _ <- checkTypeE e2
  go op t1 t2 >>= (`pureInfo` False)
  where
    go :: RelOp -> EnvType -> EnvType -> TM EnvType
    go EQU {} t1 t2 = checkSameType t1 t2
    go NE {} t1 t2 = checkSameType t1 t2
    go _ EnvInt EnvInt = pure EnvBoolean
    -- go _ EnvString EnvString = pure EnvBoolean
    go _ t1 t2 = throwError $ InvalidRelOp (hasPosition op) (toLangString t1) (toLangString t2)
    throwOfVoid :: EnvType -> EnvType -> TM EnvType
    throwOfVoid t1 t2 = throwError $ InvalidRelTypeVoid (hasPosition op) (toLangString t1) (toLangString t2)
    checkSameType :: EnvType -> EnvType -> TM EnvType
    checkSameType t1 t2 = do
      when (t1 == EnvVoid || t2 == EnvVoid) (void $ throwOfVoid t1 t2)
      c1 <- t1 `isOfType` t2
      c2 <- t2 `isOfType` t1
      if c1 || c2
        then pure EnvBoolean
        else throwError $ InvalidRelType (hasPosition op) (toLangString t1) (toLangString t2)
checkTypeE (EAnd pos e1 e2) = do
  ValueInfo t1 _ <- checkTypeE e1
  ValueInfo t2 _ <- checkTypeE e2
  if t1 == EnvBoolean && t2 == EnvBoolean
    then pureInfo EnvBoolean False
    else throwError $ InvalidAndType pos (toLangString t1) (toLangString t2)
checkTypeE (EOr pos e1 e2) = do
  ValueInfo t1 _ <- checkTypeE e1
  ValueInfo t2 _ <- checkTypeE e2
  if t1 == EnvBoolean && t2 == EnvBoolean
    then pureInfo EnvBoolean False
    else throwError $ InvalidOrType pos (toLangString t1) (toLangString t2)

validateFunctionCall :: String -> Pos -> [Expr] -> TM ValueInfo
validateFunctionCall functionName pos argsExps = do
  functionObject <- getEnvFunction functionName
  case functionObject of
    Nothing -> throwError $ UnknownFunction functionName pos
    Just EnvFunction {envFunRetType = retType, envFunArgs = args} -> do
      validateFunctionArgs functionName pos args argsExps
      pureInfo retType False

isOfType :: EnvType -> EnvType -> TM Bool
isOfType t1 t2 | t1 == t2 = pure True
isOfType t1@EnvClass {} t2@EnvClass {} = isSubtype t1 t2
isOfType _ _ = pure False

isSubtype :: EnvType -> EnvType -> TM Bool
isSubtype (EnvClass subClass) (EnvClass superClass) = do
  maybeSubClass <- getEnvClass subClass
  case maybeSubClass of
    Nothing -> throwError $ UnknownClass subClass Nothing
    Just EnvClassObject {envClassSuperClass = Nothing} -> pure False
    Just EnvClassObject {envClassSuperClass = Just oneUp} -> isOfType (EnvClass oneUp) (EnvClass superClass)
isSubtype _ _ = pure False

validateFunctionArgs :: String -> Pos -> [(EnvType, String)] -> [Expr] -> TM ()
validateFunctionArgs _ _ [] [] = pure ()
validateFunctionArgs funName pos ((expectedType, _) : expectedArgs) (arg : args) = do
  ValueInfo actualType _ <- checkTypeE arg
  ofType <- actualType `isOfType` expectedType
  if ofType
    then validateFunctionArgs funName pos expectedArgs args
    else throwError $ InvalidFunctionArgType pos funName (toLangString expectedType) (toLangString actualType)
validateFunctionArgs name pos _ _ = throwError $ InvalidFunctionArgCount pos name

getClassField :: Pos -> Expr -> String -> TM ValueInfo
getClassField pos obj name = do
  ValueInfo objType _ <- checkTypeE obj
  case (objType, name) of
    (EnvClass className, _) -> do
      maybeClass <- getEnvClass className
      case maybeClass of
        Nothing -> throwError $ UnknownClass className pos
        Just objClass -> do
          maybeField <- findFieldInClass name (Just objClass)
          case maybeField of
            Nothing -> throwError $ UnknownClassField name className pos
            Just t -> pureInfo t True
    (EnvArray _, "length") -> pureInfo EnvInt False
    _ -> throwError $ InvalidClassFieldAccess pos (toLangString objType)

checkMethodCall :: Pos -> Expr -> String -> [Expr] -> TM ValueInfo
checkMethodCall pos obj funName argsExps = do
  ValueInfo objType _ <- checkTypeE obj
  case objType of
    EnvClass className -> do
      maybeClass <- getEnvClass className
      case maybeClass of
        Nothing -> throwError $ UnknownClass className pos
        Just classObj -> do
          maybeMethod <- findMethodInClass funName (Just classObj)
          case maybeMethod of
            Nothing -> throwError $ UnknownClassMethod funName className pos
            Just EnvFunction {envFunArgs = args, envFunRetType = retType} -> do
              validateFunctionArgs funName pos args argsExps
              pureInfo retType False
    _ -> throwError $ InvalidClassMethodAccess pos (toLangString objType)

checkArrAccess :: Pos -> Expr -> Expr -> TM ValueInfo
checkArrAccess pos lVal e = do
  ValueInfo lValType _ <- checkTypeE lVal
  case lValType of
    EnvArray t -> do
      ValueInfo expType _ <- checkTypeE e
      if expType == EnvInt
        then pureInfo t True
        else throwError $ InvalidArrayIndex pos (toLangString expType)
    _ -> throwError $ InvalidArrayAccess pos (toLangString lValType)
