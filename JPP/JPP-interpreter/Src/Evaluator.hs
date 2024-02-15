module Src.Evaluator where

import Src.Jabba.Abs

import Src.Types ( VarRef (..) )
import Src.Errors
import Control.Monad.Except (ExceptT, runExceptT, throwError, void, MonadIO (liftIO))
import Control.Monad.Reader (ReaderT, runReaderT, ask, local, asks)
import Control.Monad.State (StateT, evalStateT, get, put, modify, gets)
import qualified Data.Map as Map (Map, empty, insert, lookup, fromList)
import Data.Maybe (isJust, fromMaybe)
import Text.Read (readMaybe)

-- =============================================================================

newtype Env = Env { getEnv :: Map.Map Ident Loc } deriving (Show)

emptyEnv :: Env
emptyEnv = Env Map.empty

envGet :: Ident -> Env -> IM Loc
envGet i (Env m) = case Map.lookup i m of
    Nothing -> throwError $ RuntimeError Nothing TypeCheckerDidntCatch
    Just l -> pure l

insertEnv :: Ident -> Loc -> Env -> Env
insertEnv i l (Env m) = Env $ Map.insert i l m

insertEnvMany :: [(Ident, Loc)] -> Env -> Env
insertEnvMany [] e = e
insertEnvMany ((i, l):xs) e = insertEnvMany xs $ insertEnv i l e

type ModEnv = Env -> Env

-- =============================================================================

data FunBlock
    = NormalBlock Block
    | WriteStr
    | WriteInt
    | ToInt
    | ToString
    | ExitCode
    | AssertExpr

type FnArg = (Ident, VarRef)
data Value
    = VTInt Int
    | VTBool Bool
    | VTString String
    | VTUnit
    | VTFun [FnArg] FunBlock Env
    | VTTab Int (Map.Map Int Value)

instance Show Value where
    show (VTInt i) = show i
    show (VTBool b) = show b
    show (VTString s) = show s
    show VTUnit = "()"
    show (VTFun args _ _) = "Fn(" ++ show args ++ ")"
    show (VTTab s vs) = "[" ++ show s ++ "<>" ++ show vs ++ "]"

data LoopControlFlow = Break | Continue

type RetType = (Maybe Value, Maybe LoopControlFlow)

type Loc = Int

data Store = Store { nextFree :: Loc, memory :: Map.Map Loc Value, programStr :: String }

emptyStore :: Store
emptyStore = storeWithProgram ""

storeWithProgram :: String -> Store
storeWithProgram = Store 0 Map.empty

newlock :: Store -> (Loc, Store)
newlock (Store n m s) = (n, Store (n + 1) m s)

newlockM :: IM Loc
newlockM = do
    s <- get
    let (l, s') = newlock s
    put s'
    pure l

insertStore :: Loc -> Value -> Store -> Store
insertStore l v (Store n m s) = Store n (Map.insert l v m) s

insertStoreMany :: [(Loc, Value)] -> Store -> Store
insertStoreMany [] s = s
insertStoreMany ((l, v):xs) s = insertStoreMany xs $ insertStore l v s

storeGet :: Loc -> Store -> IM Value
storeGet l (Store _ m _) = case Map.lookup l m of
    Nothing -> throwError $ RuntimeError Nothing TypeCheckerDidntCatch
    Just v -> pure v

reserveNnonM :: Int -> Store -> ([Loc], Store)
reserveNnonM n s = foldl (\ (ls, s') _ -> let (nl, s'') = newlock s' in (nl:ls, s'')) ([], s) [1..n]

reserveN :: Int -> IM [Loc]
reserveN n = do
    store <- get
    -- first we reserve memory for all variables
    let (locs, store') = reserveNnonM n store
    put store'
    pure locs

type ModStore = Store -> Store

-- =============================================================================

getTabValueAt :: BNFC'Position -> Int -> Loc -> IM Value
getTabValueAt pos i l = do
    (VTTab s m) <- storeGet l =<< get
    case Map.lookup i m of
        Nothing -> throwError $ RuntimeError pos $ ArrayIndexOutOfBounds i s
        Just v -> pure v

setTabValueAt :: BNFC'Position -> Int -> Loc -> Value -> IM ()
setTabValueAt pos i l v = do
    (VTTab s m) <- storeGet l =<< get
    if i >= s || i < 0 then throwError $ RuntimeError pos $ ArrayIndexOutOfBounds i s
    else do
        let m' = Map.insert i v m
        modify $ insertStore l (VTTab s m')

initArray :: BNFC'Position -> Int -> Value -> IM Value
initArray pos size defVal = initArrayDefs pos size $ replicate size defVal

initArrayDefs :: BNFC'Position -> Int -> [Value] -> IM Value
initArrayDefs pos size defs = do
    if size < 0 then throwError $ RuntimeError pos $ ArrayInvalidSize size
    else do
        let m = Map.fromList $ zip [0..size-1] defs
        pure $ VTTab size m

-- =============================================================================

--                  errors     local   env  modify memory io
type IM a = ExceptT ErrHolder (ReaderT Env (StateT Store IO)) a


evaluate :: Program -> String -> IO (Err ())
evaluate p s = do 
    let (env, store) = makeStdLib emptyEnv $ storeWithProgram s
    let errorsT = runExceptT (evalP p)
    let envT = runReaderT errorsT env
    evalStateT envT store



evalP :: Program -> IM ()
evalP (PProgram _ is) = void $ evalIs is



evalIs :: [Instr] -> IM RetType
evalIs [] = pure (Nothing, Nothing)
evalIs ((IExpr _ e):is) = evalE e >> evalIs is
evalIs ((IDecl _ d):is) = do
    modEnv <- evalD d
    local modEnv $ evalIs is
evalIs ((IUnit _):is) = evalIs is
evalIs ((IIncr _ v):is) = modIntVar v (+1) >> evalIs is
evalIs ((IDecr _ v):is) = modIntVar v (subtract 1) >> evalIs is
evalIs ((IAss _ v e):is) = do 
    n <- evalE e
    loc <- envGet v =<< ask
    modify $ insertStore loc n
    evalIs is
evalIs ((IRet _ e):_) = do
    v <- evalE e
    pure (Just v, Nothing)
evalIs ((IRetUnit _):_) = pure (Just VTUnit, Nothing)
evalIs ((IBreak _):_) = pure (Nothing, Just Break)
evalIs ((ICont _):_) = pure (Nothing, Just Continue)
-- evalIs ((IIfElifElse _ eb b1 [] b2):is) = evalIf eb b1 b2 is
-- evalIs ((IIfElifElse pos eb b1 ((ElseIf _ eb2 b2):xs) b3):is) = evalIf eb b1 (IBlock pos [IIfElifElse pos eb2 b2 xs b3]) is
-- evalIs ((IIfElif pos eb b1 elifs):is) = evalIs (IIfElifElse pos eb b1 elifs (IBlock pos []):is)
evalIs ((IIf _ ifStmt):is) = evalIfStmt ifStmt is
evalIs whileLoop@((IWhile _ e b):is) = do
    (VTBool cond) <- evalE e
    if cond then do
        (ret, flow) <- evalB b
        if isJust ret then pure (ret, Nothing)
        else case flow of
            Nothing -> evalIs whileLoop
            Just Break -> pure (Nothing, Nothing)
            Just Continue -> evalIs whileLoop
    else evalIs is
evalIs whileLoop@((IWhileFin _ e b bFin):is) = do
    (VTBool cond) <- evalE e
    if cond then do
        (ret, flow) <- evalB b
        if isJust ret then pure (ret, Nothing)
        else case flow of
            Nothing -> evalIs whileLoop
            Just Break -> pure (Nothing, Nothing)
            Just Continue -> evalIs whileLoop
    else evalB bFin >> evalIs is
evalIs ((IFor _ v e1 e2 b):is) = do
    (VTInt n1) <- evalE e1
    (VTInt n2) <- evalE e2
    let range = [n1..n2]
    loc <- newlockM
    local (insertEnv v loc) (runFor loc range b) >>= handleRetType is
evalIs ((IBBlock _ b):is) = evalB b >>= handleRetType is
evalIs ((DFun _ fN args _ b):is) = do
    let args' = map resolveDeclArg args
    env <- ask
    loc <- newlockM
    let f = VTFun args' (NormalBlock b) (insertEnv fN loc env)
    modify $ insertStore loc f
    local (insertEnv fN loc) (evalIs is)
evalIs ((DFunUnit pos fN args b):is) = evalIs (DFun pos fN args (TVoid pos) b:is)
evalIs ((ITabAss pos aN i e):is) = do
    n <- evalE e
    (VTInt aId) <- evalE i
    l <- envGet aN =<< ask
    setTabValueAt pos aId l n
    evalIs is



evalIfStmt :: IfStmt -> [Instr] -> IM RetType
evalIfStmt (IfIf pos e b) is = evalIf e b (IBlock pos []) is
evalIfStmt (IfElse _ e b1 b2) is = evalIf e b1 b2 is
evalIfStmt (IfElseIf pos e b1 elif) is = evalIf e b1 (IBlock pos [IIf pos elif]) is


evalIf :: Expr -> Block -> Block -> [Instr] -> IM RetType
evalIf eb b1 b2 is = do
    retTypeB <- evalE eb
    case retTypeB of
        VTBool b -> 
            if b then evalB b1 >>= handleRetType is
                else evalB b2 >>= handleRetType is
        _ -> throwError $ RuntimeError Nothing TypeCheckerDidntCatch



resolveDeclArg :: Arg -> (Ident, VarRef)
resolveDeclArg (RefMutArg    _ i _) = (i, VRRef)
resolveDeclArg (RefConstArg  _ i _) = (i, VRRef)
resolveDeclArg (CopyMutArg   _ i _) = (i, VRCopy)
resolveDeclArg (CopyConstArg _ i _) = (i, VRCopy)



runFor :: Loc -> [Int] -> Block -> IM RetType
runFor _ [] _ = pure (Nothing, Nothing)
runFor loc (n:ns) b = do
    modify $ insertStore loc (VTInt n)
    (ret, flow) <- evalB b
    if isJust ret then pure (ret, Nothing)
    else case flow of
        Nothing -> runFor loc ns b
        Just Break -> pure (Nothing, Nothing)
        Just Continue -> runFor loc ns b



handleRetType :: [Instr] -> RetType -> IM RetType
handleRetType is (Nothing, Nothing) = evalIs is
handleRetType _ (Just _, Just _) = undefined -- There should never be a return and a break/continue
handleRetType _ (ret, Nothing) = pure (ret, Nothing)
handleRetType _ (Nothing, controlFlow) = pure (Nothing, controlFlow)



modIntVar :: Ident -> (Int -> Int) -> IM ()
modIntVar v f = do
    loc <- envGet v =<< ask
    store <- get
    (VTInt n) <- storeGet loc store
    modify $ insertStore loc (VTInt $ f n)



evalD :: Decl -> IM ModEnv
evalD (DVar _ vars) = addItemsToEnv vars
evalD (DVal _ vars) = addItemsToEnv vars



addItemsToEnv :: [Item] -> IM ModEnv
addItemsToEnv vars = do
    items <- mapM resolveItem vars
    addVarsToEnv items



addVarsToEnv :: [(Ident, Value)] -> IM ModEnv
addVarsToEnv items = do
    locs <- reserveN $ length items
    -- then we insert values into memory
    mapM_ (\ ((_, v), l) -> modify $ insertStore l v) $ zip items locs
    let ids = map fst items
    -- we need to tell how to modify environment
    let envMod (Env m) = Env $ foldl (\ m' (i, l) -> Map.insert i l m') m $ zip ids locs
    pure envMod



resolveItem :: Item -> IM (Ident, Value)
resolveItem (DItemVal _ i _ e) = do
    v <- evalE e
    pure (i, v)
resolveItem (DItemAuto _ i e) = do
    v <- evalE e
    pure (i, v)
resolveItem (DItem _ i t) = pure (i, defaultValueForType t)



defaultValueForType :: Type -> Value
defaultValueForType (TInt _) = VTInt 0
defaultValueForType (TBool _) = VTBool False
defaultValueForType (TString _) = VTString ""
defaultValueForType (TVoid _) = VTUnit
defaultValueForType TFun {} = undefined -- type checker should catch this
defaultValueForType TTab {} = undefined -- type checker should catch this



evalB :: Block -> IM RetType
evalB (IBlock _ is) = evalIs is



doDivLikeOp :: (Int -> Int -> Int) -> BNFC'Position -> Int -> Int -> IM Value
doDivLikeOp op pos n1 n2 = if n2 == 0 then throwError $ RuntimeError pos ZeroDiv
    else pure $ VTInt $ op n1 n2



evalE :: Expr -> IM Value

evalE (EVarName pos v) = getVarValue pos v
evalE (EIntLit _ n) = pure $ VTInt $ fromInteger n
evalE (EBoolLitTrue _) = pure $ VTBool True
evalE (EBoolLitFalse _) = pure $ VTBool False
evalE (EUnitLiteral _) = pure VTUnit
evalE (EStringLit _ s) = pure $ VTString s
evalE (ENeg _ (ONeg _) e) = do
    (VTInt n) <- evalE e
    pure (VTInt $ negate n)
evalE (ENot _ (ONot _) e) = do
    (VTBool b) <- evalE e
    pure (VTBool $ not b)
evalE (EMul _ e1 op e2) = do
    (VTInt n1) <- evalE e1
    (VTInt n2) <- evalE e2
    case op of
        OMul _ -> pure (VTInt $ n1 * n2)
        ODiv pos -> doDivLikeOp div pos n1 n2
        OMod pos -> doDivLikeOp mod pos n1 n2
evalE (ESum pos e1 op e2) = do
    v1 <- evalE e1
    v2 <- evalE e2
    case (op, v1, v2) of
        (OPlus  _, VTInt n1, VTInt n2) -> pure (VTInt $ n1 + n2)
        (OMinus _, VTInt n1, VTInt n2) -> pure (VTInt $ n1 - n2)
        (OPlus  _, VTString s1, VTString s2) -> pure (VTString $ s1 ++ s2)
        _ -> throwError $ RuntimeError pos TypeCheckerDidntCatch
evalE (ERel _ e1 op e2) = do
    v1 <- evalE e1
    v2 <- evalE e2
    b <- doRelOp op v1 v2
    pure $ VTBool b
evalE (EBAnd _ e1 (OAnd _) e2) = do
    (VTBool b1) <- evalE e1
    if b1
        then do
            (VTBool b2) <- evalE e2
            pure $ VTBool b2
        else pure $ VTBool False
evalE (EBOr _ e1 (OOr _) e2) = do
    (VTBool b1) <- evalE e1
    if b1
        then pure $ VTBool True
        else do
            (VTBool b2) <- evalE e2
            pure $ VTBool b2
evalE (ETer _ eb e1 e2) = do
    (VTBool b) <- evalE eb
    if b
        then evalE e1
        else evalE e2
evalE (ERun pos e argsE) = do
    (VTFun args b env) <- evalE e
    argsVMl <- mapM evalERef argsE
    let (argsV, argsMl) = unzip argsVMl
    let (argsI, argsR) = unzip args
    argsL <- mapM insertArg $ zip3 argsR argsV argsMl
    let envInserter = insertEnvMany $ zip argsI argsL
    (ret, Nothing) <- local (\_ -> envInserter env) $ evalFunB pos b
    case ret of 
        Nothing -> pure VTUnit
        Just v  -> pure v
    where
        insertArg :: (VarRef, Value, Maybe Loc) -> IM Loc
        insertArg (r, v, ml) = do
            loc <- case r of
              VRRef -> do
                let Just l = ml
                pure l
              VRCopy -> do
                [l] <- reserveN 1
                pure l
            modify $ insertStore loc v
            pure loc
evalE (ELambda _ args b) = do
    let args' = map resolveDeclArg args
    asks $ VTFun args' (NormalBlock b)
evalE (ELambdaEmpty pos b) = evalE (ELambda pos [] b)
evalE (ELambdaExpr pos args e) = evalE (ELambda pos args (IBlock pos [IRet pos e]))
evalE (ELambdaEmptEpr pos e) = evalE (ELambda pos [] (IBlock pos [IRet pos e]))
evalE (ITabAcc pos v e) = do
    arrayLoc <- getVarLoc pos v
    (VTInt i) <- evalE e
    getTabValueAt pos i arrayLoc
evalE (ITabInit pos eS eD) = do
    (VTInt arrSize) <- evalE eS
    defVal <- evalE eD
    initArray pos arrSize defVal
evalE (ITabInitEls pos els) = do
    elsVals <- mapM evalE els
    initArrayDefs pos (length elsVals) elsVals



evalERef :: Expr -> IM (Value, Maybe Loc)
evalERef (EVarName _ v) = do
    loc <- envGet v =<< ask
    store <- get
    val <- storeGet loc store
    pure (val, Just loc)
evalERef e = do
    v <- evalE e
    pure (v, Nothing)



doRelOp :: RelOp -> Value -> Value -> IM Bool
doRelOp (REq _)  (VTInt n1)    (VTInt n2)    = pure $ n1 == n2
doRelOp (RNeq _) (VTInt n1)    (VTInt n2)    = pure $ n1 /= n2
doRelOp (RLt _)  (VTInt n1)    (VTInt n2)    = pure $ n1 <  n2
doRelOp (RGt _)  (VTInt n1)    (VTInt n2)    = pure $ n1 >  n2
doRelOp (RLeq _) (VTInt n1)    (VTInt n2)    = pure $ n1 <= n2
doRelOp (RGeq _) (VTInt n1)    (VTInt n2)    = pure $ n1 >= n2
doRelOp (REq _)  (VTBool b1)   (VTBool b2)   = pure $ b1 == b2
doRelOp (RNeq _) (VTBool b1)   (VTBool b2)   = pure $ b1 /= b2
doRelOp (REq _)  (VTString s1) (VTString s2) = pure $ s1 == s2
doRelOp (RNeq _) (VTString s1) (VTString s2) = pure $ s1 /= s2
doRelOp _ _ _ = throwError $ RuntimeError Nothing TypeCheckerDidntCatch



-- =============================================================================



getVarLoc :: BNFC'Position -> Ident -> IM Loc
getVarLoc _ v = do
    env <- ask
    envGet v env


getVarValue :: BNFC'Position -> Ident -> IM Value
getVarValue pos v = do
    loc <- getVarLoc pos v
    store <- get
    storeGet loc store



-- =============================================================================


evalFunB :: BNFC'Position -> FunBlock -> IM RetType
evalFunB _ (NormalBlock b) = evalB b
evalFunB pos WriteStr = do
    (VTString s) <- evalE (EVarName pos (Ident "s"))
    liftIO $ putStr s
    pure (Nothing, Nothing)
evalFunB pos WriteInt = do
    (VTInt n) <- evalE (EVarName pos (Ident "n"))
    liftIO $ putStr $ show n
    pure (Nothing, Nothing)
evalFunB pos ToString = do
    (VTInt n) <- evalE (EVarName pos (Ident "n"))
    pure (Just $ VTString $ show n, Nothing)
evalFunB pos ToInt = do
    (VTString s) <- evalE (EVarName pos (Ident "s"))
    case readMaybe s of
        Just n -> pure (Just $ VTInt n, Nothing)
        Nothing -> throwError $ RuntimeError pos $ CannotCast s "Integer"
evalFunB pos ExitCode = do
    (VTInt n) <- evalE (EVarName pos (Ident "n"))
    throwError $ ControlledExit n
evalFunB pos AssertExpr = do
    (VTBool b) <- evalE (EVarName pos (Ident "e"))
    if b
        then pure (Nothing, Nothing)
        else do
            codeFragment <- getCodeFragmentFromPos pos
            throwError $ RuntimeError pos $ AssertionFailed codeFragment


getCodeFragmentFromPos :: BNFC'Position -> IM String
getCodeFragmentFromPos (Just (line, col)) = do
    code' <- gets programStr
    let code = lines code'
    let line' = line - 1
    let col' = col - 1 + length "assert("
    let lineStr = code !! line'
    let lineBefore = fromMaybe "" $ code `getElem` (line' - 1)
    let lineAfter = fromMaybe "" $ code `getElem` (line' + 1)
    getAssertMessage line col' lineStr lineBefore lineAfter
getCodeFragmentFromPos Nothing = pure "no position info"


getElem :: [a] -> Int -> Maybe a
getElem [] _ = Nothing
getElem _ n | n < 0 = Nothing
getElem (x:_) 0 = Just x
getElem (_:xs) n = getElem xs (n - 1)

-- for string `assert(false)` and numer 6 generates message like
-- ```
-- assert(false);
--       ^
-- ```
getAssertMessage :: Int -> Int -> String -> String -> String -> IM String
getAssertMessage line col code codeBefore codeAfter = do
    let beforeCode = show line ++ ": "
    let beforeCodLen = length beforeCode
    let beforeBeforeCode = show (line - 1) ++ ": "
    let beforeBeforeCodLen = length beforeBeforeCode
    let beforeAfterCode = show (line + 1) ++ ": "
    let beforeAfterCodLen = length beforeAfterCode
    let codeBeforeMore = if line > 1 then beforeBeforeCode ++ codeBefore ++ "\n" else ""
    let codeAfterMore = beforeAfterCode ++ codeAfter
    let maxLen = maximum [beforeCodLen, beforeBeforeCodLen, beforeAfterCodLen]
    pure $ codeBeforeMore ++ beforeCode ++ code ++ "\n" ++ duplicate " " (col + maxLen) ++ "^\n" ++ codeAfterMore

duplicate :: String -> Int -> String
duplicate string n = concat $ replicate n string

writeStrDecl :: (Ident, Value)
writeStrDecl = (Ident "writeStr", VTFun [(Ident "s", VRCopy)] WriteStr emptyEnv)

writeIntDecl :: (Ident, Value)
writeIntDecl = (Ident "writeInt", VTFun [(Ident "n", VRCopy)] WriteInt emptyEnv)

toStringDecl :: (Ident, Value)
toStringDecl = (Ident "toString", VTFun [(Ident "n", VRCopy)] ToString emptyEnv)

toIntDecl :: (Ident, Value)
toIntDecl = (Ident "toInt", VTFun [(Ident "s", VRCopy)] ToInt emptyEnv)

exitCode :: (Ident, Value)
exitCode = (Ident "exit", VTFun [(Ident "n", VRCopy)] ExitCode emptyEnv)

assertExpr :: (Ident, Value)
assertExpr = (Ident "assert", VTFun [(Ident "e", VRCopy)] AssertExpr emptyEnv)


stdLib :: [(Ident, Value)]
stdLib = [writeStrDecl, writeIntDecl, toStringDecl, toIntDecl, exitCode, assertExpr]


makeStdLib :: Env -> Store -> (Env, Store)
makeStdLib env store = 
    let stdLibCount = length stdLib in
    let (locs, store') = reserveNnonM stdLibCount store in
    let (names, values) = unzip stdLib in
    let env' = insertEnvMany (zip names locs) env in
    let store'' = insertStoreMany (zip locs values) store' in
    (env', store'')
