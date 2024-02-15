module Src.TypeChecker where

import Src.Jabba.Abs
import qualified Data.Map as Map
import qualified Data.Set as Set
import Control.Monad.Except (ExceptT, MonadError (throwError), runExceptT, when)
import Control.Monad.State (State, evalState, MonadState (get, put), modify)
import Src.Errors
    ( ErrType (..),
      ErrHolder (TypeChecker),
      Err )
import Src.Types

data RetType'
    = None
    | Definitive VarType
    | Branching VarType

--             return     loop flow control
type RetType = (RetType', Bool)

newtype Env = Env (Map.Map Ident (VarType, VarMutability))

type IM a = ExceptT ErrHolder (State Env) a


emptyEnv :: Env
emptyEnv = Env Map.empty


-- Functions defined using `fun` keyword cannot be replaced using lambda expressions
-- e.g.
-- fun f() { return 1; }
-- f = || { return 2; } // Error. `fun` is equivalent to `val`
fMut :: VarMutability
fMut = VMConst


envUnion :: Env -> Env -> Env
envUnion (Env e1) (Env e2) = Env $ Map.union e1 e2


envInsert :: Ident -> (VarType, VarMutability) -> Env -> Env
envInsert i v (Env e) = Env $ Map.insert i v e


stdLib :: Env
stdLib = Env (
        Map.fromList [
            (Ident "writeStr", (Fn [(VTString, VMConst, VRCopy)] VTVoid, fMut)),
            (Ident "writeInt", (Fn [(VTInt, VMConst, VRCopy)] VTVoid, fMut)),
            (Ident "toString", (Fn [(VTInt, VMConst, VRCopy)] VTString, fMut)),
            (Ident "toInt", (Fn [(VTString, VMConst, VRCopy)] VTInt, fMut)),
            (Ident "exit", (Fn [(VTInt, VMConst, VRCopy)] VTVoid, fMut)),
            (Ident "assert", (Fn [(VTBool, VMConst, VRCopy)] VTVoid, fMut))
        ]
    )


localState :: (Env -> Env) -> IM a -> IM a
localState f m = do
    s <- get
    put (f s)
    r <- m
    put s
    return r


typeCheck :: Program -> Err ()
typeCheck = typeCheckWithEnv emptyEnv


typeCheckWithEnv :: Env -> Program -> Either ErrHolder ()
typeCheckWithEnv env p = evalState (runExceptT (checkTypeP p)) $ envUnion env stdLib


checkTypeP :: Program -> IM ()
checkTypeP (PProgram pos is) = do
    retType <- checkTypeB (IBlock pos is)
    case retType of
        (None, False) -> pure ()
        (Definitive t, _) -> throwError $ TypeChecker pos $ TopLevelProgramReturn t
        (Branching t, _) -> throwError $ TypeChecker pos $ TopLevelProgramMaybeReturn t
        (None, True) -> throwError $ TypeChecker pos TopLevelProgramLoopFlow


checkTypeI :: Instr -> IM RetType
checkTypeI (IUnit _) = pure (None, False)
checkTypeI (IExpr _ e) = checkTypeE e >> pure (None, False)
checkTypeI (IIncr pos v) = opOnVarType pos v [VTInt]
checkTypeI (IDecr pos v) = opOnVarType pos v [VTInt]
checkTypeI (IDecl _ d) = checkTypeD d
checkTypeI (IAss pos v e) = do
    (t, m) <- getVarType pos v
    et <- checkTypeE e
    if m == VMConst
        then throwError $ TypeChecker pos $ ConstantAssign v
    else if t /= et
        then throwError $ TypeChecker pos $ WrongType v t [et]
    else pure (None, False)
checkTypeI (IRet _ e) = checkTypeE e >>= (\x -> pure (x, False)) . Definitive
checkTypeI (IRetUnit _) = pure (Definitive VTVoid, False)
checkTypeI (IBreak _) = pure (None, True)
checkTypeI (ICont _) = pure (None, True)
checkTypeI (IIf _ ifStmt) = checkTypeIf ifStmt
checkTypeI (IWhile pos eb bl) = checkTypeIf (IfElse pos eb bl (IBlock pos [])) >>= \(x, _) -> pure (x, False)
checkTypeI (IWhileFin pos eb (IBlock _ bIs) (IBlock _ finIs)) = do
    retType <- checkTypeI (IWhile pos eb $ IBlock pos bIs)
    finRetType <- checkTypeB (IBlock pos finIs)
    mergeRetTypes pos retType finRetType
checkTypeI (IFor pos v e1 e2 bl) = do
    n1 <- checkTypeE e1
    n2 <- checkTypeE e2
    if n1 /= VTInt || n2 /= VTInt
        then throwError $ TypeChecker pos $ ForRangeTypeMismatch n1 n2
        else localState (envInsert v (VTInt, VMConst)) $ checkTypeI (IWhile pos (EBoolLitTrue pos) bl)
        >>= \(x, _) -> pure (x, False)
checkTypeI (IBBlock _ b) = checkTypeB b
checkTypeI (DFunUnit pos fName args b) = checkTypeI (DFun pos fName args (TVoid pos) b)
checkTypeI (DFun pos fName args t b) = do
    args' <- mapM resolveFnArg args
    checkAllNamesAreUnique Set.empty $ map (\(i, t', _, _, p) -> (i, t', p)) args'
    let fnArgs = map (\(_, t', m, r, _) -> (t', m, r)) args'
    let fnRet = absTypeToVarType t
    saveEnv <- get
    let f = Fn fnArgs fnRet
    modify $ modifyEnvForFunction [(fName, (f, fMut))] args'
    (retType, loopFlow) <- checkTypeB b
    if loopFlow then throwError $ TypeChecker pos $ FunctionLoopFlow fName
    else case retType of
        None -> when (fnRet /= VTVoid) $ throwError $ TypeChecker pos $ FunctionNoReturn fName fnRet
        Definitive t' -> when (t' /= fnRet) $ throwError $ TypeChecker pos $ FunctionReturnMismatch fName t' fnRet
        Branching t' -> if t' /= fnRet
            then throwError $ TypeChecker pos $ FunctionReturnMismatch fName t' fnRet
            else when (fnRet /= VTVoid) $ throwError $ TypeChecker pos $ FunctionMaybeReturn fName fnRet
    put saveEnv
    modify $ envUnion $ Env $ Map.fromList [(fName, (f, fMut))]
    pure (None, False)
checkTypeI (ITabAss pos v ei ev) = do
    ei' <- checkTypeE ei
    ev' <- checkTypeE ev
    if ei' /= VTInt then throwError $ TypeChecker pos $ InvalidArrayIndex ei'
    else do
        (varType, _) <- getVarType pos v
        if varType /= VTTab ev' then throwError $ TypeChecker pos $ NotAnArray v
        else pure (None, False)



checkTypeIf ::IfStmt -> IM RetType 
checkTypeIf (IfElse pos e b1 b2) = checkIf pos e b1 b2
checkTypeIf (IfIf pos e b) = checkIf pos e b (IBlock pos [])
checkTypeIf (IfElseIf pos e b nextIf) = checkIf pos e b (IBlock pos [IIf pos nextIf])



checkIf :: BNFC'Position -> Expr -> Block -> Block -> IM RetType
checkIf pos eb bl1 bl2 = do
    eb' <- checkTypeE eb
    if eb' /= VTBool
        then throwError $ TypeChecker pos $ WrongTypeOp "if else condition" eb'
        else do
            (bl1Ret, loopFlow1) <- checkTypeB bl1
            (bl2Ret, loopFlow2) <- checkTypeB bl2
            let loopFlow = loopFlow1 || loopFlow2
            case (bl1Ret, bl2Ret) of
              (None, None) -> pure None
              (Definitive vt, None) -> pure $ Branching vt
              (None, Definitive vt) -> pure $ Branching vt
              (Branching vt, None)  -> pure $ Branching vt
              (None, Branching vt)  -> pure $ Branching vt
              (Definitive vt1, Definitive vt2) -> go vt1 vt2 Definitive
              (Branching  vt1, Branching  vt2) -> go vt1 vt2 Branching
              (Definitive vt1, Branching  vt2) -> go vt1 vt2 Branching
              (Branching  vt1, Definitive vt2) -> go vt1 vt2 Branching
              >>= \x -> pure (x, loopFlow)
            where
                go :: VarType -> VarType -> (VarType -> RetType') -> IM RetType'
                go vt1 vt2 f =
                    if vt1 == vt2
                        then pure $ f vt1
                        else throwError $ TypeChecker pos $ MismatchedReturnTypes vt1 vt2


modifyEnvForFunction :: [(Ident, (VarType, VarMutability))] -> [(Ident, VarType, VarMutability, VarRef, BNFC'Position)] -> Env -> Env
modifyEnvForFunction f args = envUnion newEnv
    where
        fEnv = Env $ Map.fromList f
        vEnv = Env $ Map.fromList $ map (\(i, t, m, _, _) -> (i, (t, m))) args
        newEnv = envUnion fEnv vEnv


resolveFnArg :: Arg -> IM (Ident, VarType, VarMutability, VarRef, BNFC'Position)
resolveFnArg (RefMutArg    p v t) = pure (v, absTypeToVarType t, VMMut,   VRRef,  p)
resolveFnArg (RefConstArg  p v t) = pure (v, absTypeToVarType t, VMConst, VRRef,  p)
resolveFnArg (CopyMutArg   p v t) = pure (v, absTypeToVarType t, VMMut,   VRCopy, p)
resolveFnArg (CopyConstArg p v t) = pure (v, absTypeToVarType t, VMConst, VRCopy, p)


checkTypeB :: Block -> IM RetType
checkTypeB (IBlock _ is) = localState id $ checkTypeBHelper (None, False) is

--                  previous   next  
checkTypeBHelper :: RetType -> [Instr] -> IM RetType
checkTypeBHelper prev [] = pure prev
checkTypeBHelper prev (i:is) = do
    ir <- checkTypeI i
    let pos = hasPosition i
    merged <- mergeRetTypes pos prev ir
    checkTypeBHelper merged is


mergeRetTypes :: BNFC'Position -> RetType -> RetType -> IM RetType
mergeRetTypes pos (r1, l1) (r2, l2) = do
    let loopFlow = l1 || l2
    case (r1, r2) of
        (None, _) -> pure r2
        (_, None) -> pure r1
        (Branching  t1, Definitive t2) -> merge t1 t2 Definitive
        (Branching  t1, Branching  t2) -> merge t1 t2 Branching
        (Definitive t1, Branching  t2) -> merge t1 t2 Definitive
        (Definitive t1, Definitive t2) -> merge t1 t2 Definitive
        >>= \x -> pure (x, loopFlow)
    where
        merge :: VarType -> VarType -> (VarType -> RetType') -> IM RetType'
        merge t1 t2 f =
            if t1 == t2
                then pure $ f t1
                else throwError $ TypeChecker pos $ MismatchedReturnTypes t1 t2


checkTypeD :: Decl -> IM RetType
checkTypeD (DVar pos its) = declareNewVariables pos its VMMut
checkTypeD (DVal pos its) = declareNewVariables pos its VMConst


checkTypeItem :: Item -> IM (Ident, VarType, BNFC'Position)
checkTypeItem (DItemVal pos v t e) = do
    let t' = absTypeToVarType t
    et' <- checkTypeE e
    if t' == et'
        then pure (v, t', pos)
        else throwError $ TypeChecker pos $ WrongType v t' [et']
checkTypeItem (DItem pos f TFun {}) = throwError $ TypeChecker pos $ FunctionWithoutInitializer f
checkTypeItem (DItem pos a TTab {}) = throwError $ TypeChecker pos $ TableWithoutInitializer a
checkTypeItem (DItem pos v t) = pure (v, absTypeToVarType t, pos)
checkTypeItem (DItemAuto pos v e) = do
    et' <- checkTypeE e
    pure (v, et', pos)


checkTypeE :: Expr -> IM VarType
checkTypeE (EVarName pos v) = fst <$> getVarType pos v
checkTypeE (EIntLit _ _) = pure VTInt
checkTypeE (EStringLit _ _) = pure VTString
checkTypeE (EBoolLitFalse _) = pure VTBool
checkTypeE (EBoolLitTrue _) = pure VTBool
checkTypeE (EUnitLiteral _) = pure VTVoid
checkTypeE (ERun pos ef es) = do
    fType <- checkTypeE ef
    (args, ret) <- decunstructFnType pos fType
    callArgs <- mapM (getCallExprType pos) es
    if length args /= length callArgs
        then throwError $ TypeChecker pos $ WrongNumberOfArgs (length args) (length callArgs)
        else do
            let args' = zip3 callArgs args [0..]
            mapM_ (checkCallArgument pos) args'
            pure ret
checkTypeE (ENeg pos (ONeg _) e) = checkExprSingleOp "negation" pos e VTInt
checkTypeE (ENot pos (ONot _) e) = checkExprSingleOp "negation" pos e VTBool
checkTypeE (EMul pos _ (ODiv _) (EIntLit _ 0)) = throwError $ TypeChecker pos ZeroLiternalDiv
checkTypeE (EMul pos e1 op e2) = checkExprBiOp (getMulOpName op) pos e1 e2 [VTInt]
checkTypeE (ESum pos e1 op@(OPlus _) e2) = checkExprBiOp (getSumOpName op) pos e1 e2 [VTInt, VTString]
checkTypeE (ESum pos e1 op e2) = checkExprBiOp (getSumOpName op) pos e1 e2 [VTInt]
checkTypeE (ERel pos e1 op@(REq _) e2) = checkExprBiOp (getRelOpName op) pos e1 e2 [VTInt, VTString, VTBool] >> pure VTBool
checkTypeE (ERel pos e1 op@(RNeq _) e2) = checkExprBiOp (getRelOpName op) pos e1 e2 [VTInt, VTString, VTBool] >> pure VTBool
checkTypeE (ERel pos e1 op e2) = checkExprBiOp (getRelOpName op) pos e1 e2 [VTInt] >> pure VTBool
checkTypeE (EBAnd pos e1 (OAnd _) e2) = checkExprBiOp "and" pos e1 e2 [VTBool]
checkTypeE (EBOr pos e1 (OOr _) e2) = checkExprBiOp "or" pos e1 e2 [VTBool]
checkTypeE (ETer pos eb e1 e2) = do
    t1 <- checkTypeE e1
    t2 <- checkTypeE e2
    tb <- checkTypeE eb
    if tb /= VTBool
        then throwError $ TypeChecker pos $ WrongTypeOp "ternary operator condition" tb
        else if t1 /= t2
            then throwError $ TypeChecker pos $ TernaryMismatch t1 t2
            else pure t1
checkTypeE (ELambda pos args b) = do
    let fName = Ident "lambda"
    args' <- mapM resolveFnArg args
    checkAllNamesAreUnique Set.empty $ map (\(i, t, _, _, p) -> (i, t, p)) args'
    let fnArgs = map (\(_, t, m, r, _) -> (t, m, r)) args'
    saveEnv <- get
    modify $ modifyEnvForFunction [] args'
    (retType, loopFlow) <- checkTypeB b
    deducedType <- if loopFlow 
        then throwError $ TypeChecker pos $ FunctionLoopFlow fName
        else case retType of
            None -> pure VTVoid  -- no return -> void
            Definitive t -> pure t
            Branching t -> if t /= VTVoid 
                then throwError $ TypeChecker pos $ FunctionMaybeReturn fName t
                else pure t
    put saveEnv
    pure $ Fn fnArgs deducedType
checkTypeE (ELambdaExpr _ args e) = do
    args' <- mapM resolveFnArg args
    checkAllNamesAreUnique Set.empty $ map (\(i, t, _, _, p) -> (i, t, p)) args'
    let fnArgs = map (\(_, t, m, r, _) -> (t, m, r)) args'
    saveEnv <- get
    modify $ modifyEnvForFunction [] args'
    deducedType <- checkTypeE e
    put saveEnv
    pure $ Fn fnArgs deducedType
checkTypeE (ELambdaEmpty pos b) = checkTypeE (ELambda pos [] b)
checkTypeE (ELambdaEmptEpr pos b) = checkTypeE (ELambdaExpr pos [] b)
checkTypeE (ITabAcc pos v e) = do
    ei' <- checkTypeE e
    if ei' /= VTInt then throwError $ TypeChecker pos $ InvalidArrayIndex ei'
    else do
        (vType, _) <- getVarType pos v
        deconstructTabType pos vType v
checkTypeE (ITabInit pos el ev) = do
    el' <- checkTypeE el
    ev' <- checkTypeE ev
    if el' /= VTInt then throwError $ TypeChecker pos $ InvalidArrayIndex el'
    else pure $ VTTab ev'
checkTypeE (ITabInitEls pos els) = do
    els' <- mapM checkTypeE els
    if null els' then throwError $ TypeChecker pos EmptyArray
    else do
        let elType = head els'
        if any (/= elType) els' then throwError $ TypeChecker pos $ ArrayElementsMismatch els'
        else pure $ VTTab elType


getMulOpName :: MulOp -> String
getMulOpName (OMul _) = "multiplication"
getMulOpName (ODiv _) = "division"
getMulOpName (OMod _) = "modulo"

getSumOpName :: PlsOp -> String
getSumOpName (OPlus _)  = "addition"
getSumOpName (OMinus _) = "subtraction"

getRelOpName :: RelOp -> String
getRelOpName (REq _)  = "equality"
getRelOpName (RNeq _) = "inequality"
getRelOpName (RLt _)  = "less than"
getRelOpName (RGt _)  = "greater than"
getRelOpName (RLeq _) = "less or equal"
getRelOpName (RGeq _) = "greater or equal"


getVarType :: BNFC'Position -> Ident -> IM (VarType, VarMutability)
getVarType pos v = do
    Env env <- get
    case Map.lookup v env of
        Just (t, m) -> pure (t, m)
        Nothing -> throwError $ TypeChecker pos $ NotDeclVar v


getF :: BNFC'Position -> Ident -> IM ([FnArg], VarType)
getF pos fn = do
    Env env <- get
    case Map.lookup fn env of
        Just (Fn args ret, _) -> pure (args, ret)
        _ -> throwError $ TypeChecker pos $ NotDeclFun fn


decunstructFnType :: BNFC'Position -> VarType -> IM ([FnArg], VarType)
decunstructFnType _ (Fn args ret) = pure (args, ret)
decunstructFnType pos t = throwError $ TypeChecker pos $ NotFnType t


deconstructTabType :: BNFC'Position -> VarType -> Ident -> IM VarType
deconstructTabType _ (VTTab t) _ = pure t
deconstructTabType pos _ v = throwError $ TypeChecker pos $ NotAnArray v


getCallExprType :: BNFC'Position -> Expr -> IM (VarType, Maybe VarMutability)
getCallExprType _ (EVarName _ v) = do
    (t, m) <- getVarType undefined v
    pure (t, Just m)
getCallExprType _ e = do
    t <- checkTypeE e
    pure (t, Nothing)


checkCallArgument :: BNFC'Position -> ((VarType, Maybe VarMutability), (VarType, VarMutability, VarRef), Int) -> IM ()
checkCallArgument pos ((t, mm), (t', m, r), idd) = do
    checkCallArgumentType pos t t' idd
    checkCallArgumentMutability pos mm (m, r) idd


checkCallArgumentType :: BNFC'Position -> VarType -> VarType -> Int -> IM ()
checkCallArgumentType pos t t' idd = do
    if t == t'
        then pure ()
        else throwError $ TypeChecker pos $ WrongTypeArg idd t t'


-- There are 12 cases to check, but only 3 are problematic.
checkCallArgumentMutability :: BNFC'Position -> Maybe VarMutability -> (VarMutability, VarRef) -> Int -> IM ()
checkCallArgumentMutability pos Nothing (_, VRRef) idd = throwError $ TypeChecker pos $ ExprRefPass idd
checkCallArgumentMutability pos (Just VMConst) (VMMut, VRRef) idd = throwError $ TypeChecker pos $ ImmutMutPass idd
checkCallArgumentMutability _ _ _ _ = pure ()


opOnVarType :: BNFC'Position -> Ident -> [VarType] -> IM RetType
opOnVarType pos v _ = do
    (t, m) <- getVarType pos v
    case t of
        VTInt -> case m of
            VMMut -> pure (None, False)
            VMConst -> throwError $ TypeChecker pos $ ImmutVar v
        _ -> throwError $ TypeChecker pos $ WrongType v t [VTInt]


checkAllNamesAreUnique :: Set.Set Ident -> [(Ident, VarType, BNFC'Position)] -> IM ()
checkAllNamesAreUnique _ [] = pure ()
checkAllNamesAreUnique s ((v, _, pos):xs) = do
    if Set.member v s
        then throwError $ TypeChecker pos $ VarAlreadyDecl v
        else checkAllNamesAreUnique (Set.insert v s) xs


declareNewVariables :: BNFC'Position -> [Item' BNFC'Position] -> VarMutability -> IM RetType
declareNewVariables _ its m = do
    items <- mapM checkTypeItem its
    checkAllNamesAreUnique Set.empty items
    mapM_ (\(v, t, _) -> modify (envInsert v (t, m))) items
    pure (None, False)


checkExprSingleOp :: String -> BNFC'Position -> Expr' BNFC'Position -> VarType -> IM VarType
checkExprSingleOp opN pos e t = do
    t' <- checkTypeE e
    if t == t'
        then pure t
        else throwError $ TypeChecker pos $ WrongTypeOp opN t'


checkExprBiOp :: String -> BNFC'Position -> Expr' BNFC'Position -> Expr' BNFC'Position -> [VarType] -> IM VarType
checkExprBiOp opN pos e1 e2 ts = do
    t1 <- checkTypeE e1
    t2 <- checkTypeE e2
    if t1 == t2 && t1 `elem` ts
        then pure t1
        else throwError $ TypeChecker pos $ WrongTypeBiOp opN t1 t2
