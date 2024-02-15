module Backend.Codegen where

import Backend.LLVMCode (PhiReplacer (replacePhi), replace)
import qualified Backend.LLVMCode as L
import Control.Monad.Extra (concatMapM, mapMaybeM)
import Control.Monad.State
import qualified Data.Bifunctor
import Data.List (isSuffixOf)
import qualified Data.Map as Map
import qualified Data.Set as Set
import Debug.Trace (trace)
import qualified Frontend.AbstractSyntaxTree as AST
import Frontend.TypeEvaluator (defaultValueOfType)

type CGM a = StateT CGState IO a

type BlockID = String

type PhiID = String

data StringPool = StringPool
  { spStrings :: Map.Map String String,
    spNextString :: Int
  }

emptyStringPool :: StringPool
emptyStringPool = StringPool Map.empty 0

data PhiObj = PhiObj
  { phiId :: PhiID,
    phiBlock :: BlockID,
    phiOperands :: [(BlockID, String)],
    phiHeldType :: L.Type
  }
  deriving (Eq, Ord, Show)

instance PhiReplacer PhiObj where
  replacePhi (old, new) phiObj = do
    let operands = phiOperands phiObj
    let operands' = map (Data.Bifunctor.second (replace old new)) operands
    PhiObj (phiId phiObj) (phiBlock phiObj) operands' (phiHeldType phiObj)

data BlockObj = BlockObj
  { blockId :: BlockID,
    blockCode :: [L.Instructions],
    blockTerminator :: L.Instructions
  }
  deriving (Eq, Ord, Show)

instance PhiReplacer BlockObj where
  replacePhi (old, new) blockObj = do
    let code = blockCode blockObj
    let code' = map (replacePhi (old, new)) code
    let terminator = blockTerminator blockObj
    let terminator' = replacePhi (old, new) terminator
    BlockObj (blockId blockObj) code' terminator'

data LCSEOps
  = LCSEAdd
  | LCSESub
  | LCSEMul
  | LCSEDiv
  | LCSEMod
  deriving (Eq, Ord, Show)

type LCSEState = (LCSEOps, String, String)

data CGState = CGState
  { cgNextLabel :: Int,
    cgNextRegister :: Int,
    cgNextPhi :: Int,
    cgNextForIdx :: Int,
    -- generated program
    cgProgram :: L.Program,
    -- currently generated blocks
    cgGeneratedBlocks :: [BlockObj],
    -- currently generated code
    cgCurrentBlockCode :: [L.Instructions],
    cgCurrentBlock :: BlockID,
    cgCurrentBlockTerminator :: Maybe L.Instructions,
    cgCurrentLCSE :: Map.Map LCSEState String,
    -- variable definitions in blocks
    -- (variable, block) |-> operand (register or value)
    cgCurrentDef :: Map.Map String (Map.Map BlockID String),
    -- saved strings
    cgStringPool :: StringPool,
    -- incomplete phis
    cgIncompletePhi :: Map.Map BlockID (Map.Map String PhiID),
    -- sealed blocks
    cgSealedBlocks :: Set.Set BlockID,
    -- CFG reverse edges
    cgBlockPreds :: Map.Map BlockID [BlockID],
    -- information about phis
    cgPhis :: Map.Map PhiID PhiObj,
    -- information about phis in blocks
    cgBlocksPhis :: Map.Map BlockID [PhiID],
    -- return type of current function
    cgCurrentFunctionRet :: L.Type,
    -- AST env
    cgAstEnv :: AST.ASTEnv
  }

instance PhiReplacer CGState where
  replacePhi (old, new) state = do
    let program = cgProgram state
    let program' = replacePhi (old, new) program
    let generatedBlocks = cgGeneratedBlocks state
    let generatedBlocks' = map (replacePhi (old, new)) generatedBlocks
    let currentBlockCode = cgCurrentBlockCode state
    let currentBlockCode' = map (replacePhi (old, new)) currentBlockCode
    let currentBlockTerminator = cgCurrentBlockTerminator state
    let currentBlockTerminator' = fmap (replacePhi (old, new)) currentBlockTerminator
    let currentDef = cgCurrentDef state
    let currentDef' = Map.map (Map.map (replace old new)) currentDef
    let incompletePhi = cgIncompletePhi state
    let incompletePhi' = Map.map (Map.map (replace old new)) incompletePhi
    let phis = cgPhis state
    let phis' = Map.map (replacePhi (old, new)) phis
    state
      { cgProgram = program',
        cgGeneratedBlocks = generatedBlocks',
        cgCurrentBlockCode = currentBlockCode',
        cgCurrentBlockTerminator = currentBlockTerminator',
        cgCurrentDef = currentDef',
        cgIncompletePhi = incompletePhi',
        cgPhis = phis'
      }

runStatePhiReplacer :: String -> String -> CGM ()
runStatePhiReplacer old new = do
  modify $ replacePhi (old, new)

emptyCGState :: AST.ASTEnv -> CGState
emptyCGState env =
  CGState
    { cgNextLabel = 0,
      cgNextRegister = 0,
      cgNextPhi = 0,
      cgNextForIdx = 0,
      cgProgram = L.Program [] [] [],
      cgGeneratedBlocks = [],
      cgCurrentBlockCode = [],
      cgCurrentBlock = "",
      cgCurrentBlockTerminator = Nothing,
      cgCurrentDef = Map.empty,
      cgStringPool = emptyStringPool,
      cgIncompletePhi = Map.empty,
      cgSealedBlocks = Set.empty,
      cgBlockPreds = Map.empty,
      cgPhis = Map.empty,
      cgBlocksPhis = Map.empty,
      cgCurrentFunctionRet = L.Void,
      cgAstEnv = env,
      cgCurrentLCSE = Map.empty
    }

emit :: L.Instructions -> CGM ()
emit instruction = do
  currentTerminator <- gets cgCurrentBlockTerminator
  case currentTerminator of
    Nothing -> do
      state <- get
      let code = cgCurrentBlockCode state
      put $ state {cgCurrentBlockCode = code ++ [instruction]}
    Just _ -> pure ()

nextLabel :: String -> CGM String
nextLabel label = do
  state <- get
  let nextLabel = cgNextLabel state
  put $ state {cgNextLabel = nextLabel + 1}
  pure $ label ++ "." ++ show nextLabel

nextForIdx :: String -> CGM String
nextForIdx label = do
  state <- get
  let nextForIdx = cgNextForIdx state
  put $ state {cgNextForIdx = nextForIdx + 1}
  pure $ label ++ ".for.idx." ++ show nextForIdx

data RegisterType
  = IntR
  | BoolR
  | StringR
  | PtrR
  deriving (Eq, Ord, Show)

registerTypeOfType :: L.Type -> RegisterType
registerTypeOfType L.Int = IntR
registerTypeOfType L.Bool = BoolR
registerTypeOfType L.String = StringR
registerTypeOfType _ = PtrR

unwrapRegisterType :: RegisterType -> String
unwrapRegisterType IntR = "r"
unwrapRegisterType BoolR = "b"
unwrapRegisterType StringR = "s"
unwrapRegisterType PtrR = "p"

nextRegister :: RegisterType -> CGM String
nextRegister ty = do
  state <- get
  let nextRegister = cgNextRegister state
  let unwrapped = unwrapRegisterType ty
  put $ state {cgNextRegister = nextRegister + 1}
  pure $ "%" ++ unwrapped ++ show nextRegister

nextPhi :: CGM String
nextPhi = nextPhiPrefixed "phi"

nextPhiPrefixed :: String -> CGM String
nextPhiPrefixed prefix = do
  state <- get
  let nextPhi = cgNextPhi state
  put $ state {cgNextPhi = nextPhi + 1}
  pure $ "%" ++ prefix ++ show nextPhi

getPhiObj :: PhiID -> CGM PhiObj
getPhiObj phiId = do
  phis <- gets cgPhis
  trace ("getPhiObj: " ++ phiId) $ pure ()
  pure $ Map.findWithDefault undefined phiId phis

addToStringPool :: String -> CGM String
addToStringPool s = do
  sp <- gets cgStringPool
  let strings = spStrings sp
  case Map.lookup s strings of
    Just i -> pure i
    Nothing -> do
      let nextStringNumber = spNextString sp
      let nextString = "@str." ++ show nextStringNumber
      let newStrings = Map.insert s nextString strings
      let newStringPool = sp {spStrings = newStrings, spNextString = nextStringNumber + 1}
      modify $ \state -> state {cgStringPool = newStringPool}
      addStringToLiterals nextString s
      pure nextString

addStringToLiterals :: String -> String -> CGM ()
addStringToLiterals name value = do
  program <- gets cgProgram
  let strings = L.progStrings program
  let newStrings = strings ++ [L.StringLiteral name value]
  let newProgram = program {L.progStrings = newStrings}
  modify $ \state -> state {cgProgram = newProgram}

localState :: (CGState -> CGState) -> CGM a -> CGM a
localState f action = do
  state <- get
  put $ f state
  result <- action
  put state
  pure result

makePhi :: PhiObj -> CGM L.Instructions
makePhi PhiObj {phiId = phiId, phiHeldType = heldType, phiOperands = operands} = do
  trace ("makePhi: " ++ phiId) $ pure ()
  pure $ L.Phi phiId heldType operands

makeBlock :: BlockObj -> CGM [L.Instructions]
makeBlock BlockObj {blockId = blockId, blockCode = code, blockTerminator = terminator} = do
  trace ("makeBlock: " ++ blockId) $ pure ()
  currentPhis <- gets cgBlocksPhis
  trace ("makeBlock: " ++ blockId ++ " currentPhis:" ++ show currentPhis) $ pure ()
  phisIds <- gets $ (Map.! blockId) . cgBlocksPhis
  phis <- mapM getPhiObj phisIds
  phisCode <- mapM makePhi phis
  let label = L.Label blockId
  let code' = [label] ++ phisCode ++ code ++ [terminator]
  pure code'

saveFunction :: String -> [(String, L.Type)] -> L.Type -> CGM ()
saveFunction name args ret = do
  generatedBlocks <- gets cgGeneratedBlocks
  blocksCode <- concatMapM makeBlock generatedBlocks
  let function = L.Function name args ret blocksCode
  program <- gets cgProgram
  addFunctionToState function
  cleanAfterFunction

addFunctionToState :: L.Function -> CGM ()
addFunctionToState function = do
  program <- gets cgProgram
  let currentFunctions = L.progFunctions program
  let newFunctions = currentFunctions ++ [function]
  modify $ \state -> state {cgProgram = program {L.progFunctions = newFunctions}}

cleanAfterFunction :: CGM ()
cleanAfterFunction = do
  modify $ \state -> state {cgGeneratedBlocks = []}
  modify $ \state -> state {cgCurrentBlockCode = []}
  modify $ \state -> state {cgCurrentBlock = ""}
  modify $ \state -> state {cgCurrentBlockTerminator = Nothing}
  modify $ \state -> state {cgCurrentDef = Map.empty}
  modify $ \state -> state {cgIncompletePhi = Map.empty}
  modify $ \state -> state {cgSealedBlocks = Set.empty}
  modify $ \state -> state {cgBlockPreds = Map.empty}
  modify $ \state -> state {cgPhis = Map.empty}
  modify $ \state -> state {cgBlocksPhis = Map.empty}
  modify $ \state -> state {cgCurrentLCSE = Map.empty}

codegen :: AST.AST -> AST.ASTEnv -> IO L.Program
codegen ast astEnv = do
  ((), CGState {cgProgram = program}) <- runStateT (codegenAST ast) $ emptyCGState astEnv
  pure program

cleanBlockWithName :: String -> CGM ()
cleanBlockWithName name = do
  modify $ \state ->
    state
      { cgBlocksPhis = Map.insertWith (\_ x -> x) name [] $ cgBlocksPhis state
      }
  modify $ \state ->
    state
      { cgCurrentBlock = name,
        cgCurrentBlockCode = [],
        cgCurrentBlockTerminator = Nothing,
        cgCurrentLCSE = Map.empty
      }

saveLCSE :: LCSEOps -> String -> String -> String -> CGM ()
saveLCSE op left right result = do
  let state = (op, left, right)
  currentLCSE <- gets cgCurrentLCSE
  let newLCSE = Map.insert state result currentLCSE
  modify $ \state -> state {cgCurrentLCSE = newLCSE}

readLCSE :: LCSEOps -> String -> String -> CGM (Maybe String)
readLCSE op left right = do
  let state = (op, left, right)
  currentLCSE <- gets cgCurrentLCSE
  pure $ Map.lookup state currentLCSE

emitCurrentBlockAndStartNew :: String -> CGM ()
emitCurrentBlockAndStartNew newBlockLabel = do
  closeCurrentBlock False
  cleanBlockWithName newBlockLabel

closeCurrentBlock :: Bool -> CGM ()
closeCurrentBlock canBeTriviallyClosed = do
  currentBlockLabel <- gets cgCurrentBlock
  currentBlockCode <- gets cgCurrentBlockCode
  currentBlockTerminator <- gets cgCurrentBlockTerminator
  currentReturnType <- gets cgCurrentFunctionRet
  preds <- getBlockPreds currentBlockLabel
  case currentBlockTerminator of
    Nothing | canBeTriviallyClosed -> terminate currentBlockLabel currentBlockCode L.RetVoid
    Nothing | null preds -> pure ()
    Nothing -> terminate currentBlockLabel currentBlockCode $ L.Ret currentReturnType "undef"
    Just terminator -> terminate currentBlockLabel currentBlockCode terminator
  where
    terminate :: BlockID -> [L.Instructions] -> L.Instructions -> CGM ()
    terminate currentBlockLabel currentBlockCode terminator = do
      let blockObj = BlockObj currentBlockLabel currentBlockCode terminator
      modify $ \state -> state {cgGeneratedBlocks = cgGeneratedBlocks state ++ [blockObj]}

addEdgeInCFG :: String -> CGM ()
addEdgeInCFG nextBlock = do
  currentBlock <- gets cgCurrentBlock
  trace ("addEdgeInCFG: " ++ currentBlock ++ " -> " ++ nextBlock) $ pure ()
  blockPreds <- gets cgBlockPreds
  let newBlockPreds = addEdgeInGraph nextBlock currentBlock blockPreds
  modify $ \state -> state {cgBlockPreds = newBlockPreds}

addEdgeInGraph :: BlockID -> BlockID -> (Map.Map BlockID [BlockID] -> Map.Map BlockID [BlockID])
addEdgeInGraph from to currentMap = do
  let fromSuccs = Map.findWithDefault [] from currentMap
  let newFromSuccs = fromSuccs ++ [to]
  Map.insert from newFromSuccs currentMap

appendBrIfNotTerminated :: String -> CGM ()
appendBrIfNotTerminated blockLabel = do
  currentBlockTerminator <- gets cgCurrentBlockTerminator
  currentBlock <- gets cgCurrentBlock
  trace ("appendBrIfNotTerminated: " ++ blockLabel ++ " for current=" ++ show currentBlock) $ pure ()
  clearIfDeadBlockElse $ do
    case currentBlockTerminator of
      Nothing -> do
        trace ("appendBrIfNotTerminated: " ++ blockLabel ++ " NO terminator for current=" ++ show currentBlock) $ pure ()
        modify $ \state -> state {cgCurrentBlockTerminator = Just $ L.Br blockLabel}
        addEdgeInCFG blockLabel
      Just _ -> do
        trace ("appendBrIfNotTerminated: " ++ blockLabel ++ " has terminator for current=" ++ show currentBlock) $ pure ()
        pure ()

appendCbrIfNotTerminated :: String -> String -> String -> CGM ()
appendCbrIfNotTerminated cond trueBlock falseBlock = do
  trace ("appendCbrIfNotTerminated: " ++ cond ++ " " ++ trueBlock ++ " " ++ falseBlock) $ pure ()
  currentBlockTerminator <- gets cgCurrentBlockTerminator
  clearIfDeadBlockElse $ do
    case currentBlockTerminator of
      Nothing -> do
        trace ("appendCbrIfNotTerminated: " ++ cond ++ " " ++ trueBlock ++ " " ++ falseBlock ++ " no terminator") $ pure ()
        modify $ \state -> state {cgCurrentBlockTerminator = Just $ L.Cbr cond trueBlock falseBlock}
        addEdgeInCFG trueBlock
        addEdgeInCFG falseBlock
      Just _ -> do
        trace ("appendCbrIfNotTerminated: " ++ cond ++ " " ++ trueBlock ++ " " ++ falseBlock ++ " has terminator") $ pure ()
        pure ()

appendRetIfNotTerminated :: L.Type -> String -> CGM ()
appendRetIfNotTerminated ty value = do
  currentBlockTerminator <- gets cgCurrentBlockTerminator
  clearIfDeadBlockElse $ do
    case currentBlockTerminator of
      Nothing -> do
        modify $ \state -> state {cgCurrentBlockTerminator = Just $ L.Ret ty value}
      Just _ -> pure ()

isEntryBlock :: BlockID -> Bool
isEntryBlock blockId =
  let entrySuffix = ".entry"
   in entrySuffix `isSuffixOf` blockId

clearIfDeadBlockElse :: CGM () -> CGM ()
clearIfDeadBlockElse actionIfNotDead = do
  currentBlock <- gets cgCurrentBlock
  if isEntryBlock currentBlock
    then do
      trace ("clearIfDeadBlock: " ++ currentBlock ++ " is entry") $ pure ()
      actionIfNotDead
    else do
      trace ("clearIfDeadBlock: " ++ currentBlock ++ " is not entry") $ pure ()
      allPreds <- gets cgBlockPreds
      let currentBlockPreds = Map.findWithDefault [] currentBlock allPreds
      if null currentBlockPreds
        then do
          trace ("clearIfDeadBlock: " ++ currentBlock ++ " is dead") $ pure ()
          modify $ \state -> state {cgCurrentBlockCode = []}
          modify $ \state -> state {cgCurrentBlockTerminator = Nothing}
        else do
          trace ("clearIfDeadBlock: " ++ currentBlock ++ " is not dead") $ pure ()
          actionIfNotDead

-- codegen --

isClass :: AST.Def -> Bool
isClass AST.Class {} = True
isClass _ = False

codegenAST :: AST.AST -> CGM ()
codegenAST (AST.Program defs) = do
  let classes = filter isClass defs
  let functions = filter (not . isClass) defs
  mapM_ codegenDef classes
  mapM_ codegenAllMethods classes
  mapM_ codegenDef functions

codegenClass :: AST.Def -> AST.ASTClass -> CGM ()
codegenClass (AST.Class name' methods' fields') (AST.ASTClassObj name fields methods _) = do
  unless (name' == name) $ error "Class name mismatch"
  let implementors = map (\(AST.ASTMethodObj _ _ implementor) -> implementor) methods
  vMethods <- mapM codegenMakeVMethod $ zip methods implementors
  let implementations = map makeOverrideMethodName methods
  makeConstructor name fields
  madeFields <- mapM makeField fields
  let classObj = L.ClassObj name madeFields vMethods implementations
  saveClass classObj
  where
    makeOverrideMethodName :: AST.ASTMethod -> String
    makeOverrideMethodName (AST.ASTMethodObj (AST.ASTFunObj name args retType) _ implementor) =
      implementor ++ ".method." ++ takeOnlyMethodName name
    takeOnlyMethodName :: String -> String
    takeOnlyMethodName = reverse . takeWhile (/= '.') . reverse
    makeField :: AST.ASTField -> CGM (L.Type, String)
    makeField (AST.ASTFieldObj name ty) = do
      let ty' = L.fromAstType ty
      pure (ty', name)
codegenClass _ _ = error "Don't invoke"

codegenAllMethods :: AST.Def -> CGM ()
codegenAllMethods (AST.Class name methods' fields') = do
  mapM_ (codegenMakeMethod name . addSelfArgument name) methods'
  where
    addSelfArgument :: String -> AST.FunctionObj -> AST.FunctionObj
    addSelfArgument className (AST.FunctionObj name args ret block) =
      AST.FunctionObj name (AST.Argument "self" (AST.ASTClass className) : args) ret block
codegenAllMethods _ = error "Don't invoke"

iterateLikeInPython :: Monad m => (Int -> a -> m b) -> [a] -> m [b]
iterateLikeInPython f = go 1
  where
    go _ [] = pure []
    go i (x : xs) = do
      y <- f i x
      ys <- go (i + 1) xs
      pure $ y : ys

makeConstructor :: String -> [AST.ASTField] -> CGM ()
makeConstructor className fields = do
  let ret' = L.Class className
  let self = "self"
  let name = className ++ ".constructor"
  modify $ \state -> state {cgCurrentFunctionRet = ret'}
  let entryLabel = name ++ ".entry"
  cleanBlockWithName entryLabel
  sizeofTypeR <- sizeofType $ AST.ASTClass className
  emit $ L.Call "%self.i8" L.voidPtrType "_memory_allocate" [(L.Int, sizeofTypeR)]
  emit $ L.BitcastT "%self" L.voidPtrType "%self.i8" (L.Class className)
  vTableR <- nextRegister PtrR
  emit $ L.Gep vTableR (L.UnsafeRaw $ "%" ++ className) "%self" ["0", "0"]
  emit $
    L.Store
      (L.UnsafeRaw $ "%" ++ className ++ ".vtable*")
      ("@" ++ className ++ ".vtable.data")
      (L.UnsafeRaw $ "%" ++ className ++ ".vtable**")
      vTableR
  iterateLikeInPython (makeField $ "%" ++ self) fields
  appendRetIfNotTerminated ret' $ "%" ++ self
  closeCurrentBlock True
  saveFunction name [] ret'
  where
    makeField :: String -> Int -> AST.ASTField -> CGM ()
    makeField self idx (AST.ASTFieldObj name ty) = do
      let ty' = L.fromAstType ty
      let value = defaultValueOfType ty
      let classType = L.UnsafeRaw $ "%" ++ className
      defValueR <- codegenExpr value
      valueR <- nextRegister PtrR
      emit $ L.Gep valueR classType self ["0", show idx]
      emit $ L.Store ty' defValueR (L.Ptr ty') valueR
      emit $ L.Comment $ "^ field " ++ name

codegenMakeVMethod :: (AST.ASTMethod, String) -> CGM L.Type
codegenMakeVMethod (AST.ASTMethodObj (AST.ASTFunObj name args retType) _ _, className) = do
  pure $ L.Fun (L.fromAstType retType) (L.Class className : map L.fromAstType args)

codegenMakeMethod :: String -> AST.FunctionObj -> CGM ()
codegenMakeMethod className fun@(AST.FunctionObj name args retType block) = do
  codegenFunction fun

saveClass :: L.Class -> CGM ()
saveClass classObj = do
  program <- gets cgProgram
  let classes = L.progClasses program
  let newClasses = classes ++ [classObj]
  modify $ \state -> state {cgProgram = program {L.progClasses = newClasses}}

codegenDef :: AST.Def -> CGM ()
codegenDef (AST.Function funObj) = codegenFunction funObj
codegenDef classObj@AST.Class {AST.defClassName = className} = do
  classes <- gets $ AST.astEnvClasses . cgAstEnv
  let classObj' = classes Map.! className
  codegenClass classObj classObj'

codegenFunction :: AST.FunctionObj -> CGM ()
codegenFunction (AST.FunctionObj name args ret block) = do
  let ret' = L.fromAstType ret
  modify $ \state -> state {cgCurrentFunctionRet = ret'}
  let entryLabel = name ++ ".entry"
  writeVariablesToBlock entryLabel $ map (\(AST.Argument n t) -> (n, "%" ++ n)) args
  cleanBlockWithName entryLabel
  codegenBlock block
  closeCurrentBlock (ret == AST.ASTVoid)
  tryRemoveReallyTrivialBlocks
  let args' = map (\(AST.Argument n t) -> (n, L.fromAstType t)) args
  saveFunction name args' ret'

codegenBlock :: AST.Block -> CGM ()
codegenBlock (AST.Block stmts) = mapM_ codegenStmt stmts

codegenStmt :: AST.Stmt -> CGM ()
codegenStmt (AST.SBlock block) = codegenBlock block
codegenStmt (AST.SVarDecl varName ty expr@AST.Expr {AST.exprType = rValType}) = do
  expr'' <- codegenExpr expr
  expr' <-
    if ty == rValType
      then do
        pure expr''
      else do
        expr''' <- nextRegister $ registerTypeOfType $ L.fromAstType ty
        emit $ L.BitcastT expr''' (L.fromAstType rValType) expr'' (L.fromAstType ty)
        pure expr'''
  emit $ L.Comment $ "^ var decl " ++ varName
  block <- gets cgCurrentBlock
  writeVariable varName block expr'
codegenStmt (AST.SAssign lVal@AST.Expr {AST.exprType = expectedType} rVal@AST.Expr {AST.exprType = variableType}) = do
  rVal'' <- codegenExpr rVal
  trace ("codegenStmt: " ++ show lVal) $ pure ()
  rVal' <-
    if expectedType == variableType
      then do
        pure rVal''
      else do
        rVal''' <- nextRegister $ registerTypeOfType $ L.fromAstType expectedType
        emit $ L.BitcastT rVal''' (L.fromAstType variableType) rVal'' (L.fromAstType expectedType)
        pure rVal'''
  case lVal of
    AST.Expr _ (AST.EVariable name) -> do
      emit $ L.Comment $ "^ assign " ++ name
      block <- gets cgCurrentBlock
      writeVariable name block rVal'
    AST.Expr t _ -> do
      lVal' <- codegenLVal lVal
      let t' = L.fromAstType t
      let t'ptr = L.Ptr t'
      emit $ L.Comment $ "^ assign " ++ lVal'
      emit $ L.Store t' rVal' t'ptr lVal'
codegenStmt (AST.SReturn retType expr@(AST.Expr ty _)) = do
  trace (show retType) $ pure ()
  expr'' <- codegenExpr expr
  expr' <- if ty == retType then do
    pure expr''
  else do
    expr''' <- nextRegister $ registerTypeOfType $ L.fromAstType ty
    emit $ L.BitcastT expr''' (L.fromAstType ty) expr'' (L.fromAstType retType)
    pure expr'''
  appendRetIfNotTerminated (L.fromAstType retType) expr'
codegenStmt AST.SReturnVoid = do
  appendRetIfNotTerminated L.Void ""
codegenStmt (AST.SIf condExpr trueBlock falseBlock) = do
  trueLabel <- nextLabel "if.true"
  falseLabel <- nextLabel "if.false"
  endLabel <- nextLabel "if.end"

  codegenBoolExpr condExpr trueLabel falseLabel
  sealBlock trueLabel
  sealBlock falseLabel

  emitCurrentBlockAndStartNew trueLabel
  codegenBlock trueBlock
  appendBrIfNotTerminated endLabel

  emitCurrentBlockAndStartNew falseLabel
  codegenBlock falseBlock
  appendBrIfNotTerminated endLabel

  sealBlock endLabel
  emitCurrentBlockAndStartNew endLabel
codegenStmt (AST.SWhile condExpr block) = do
  condLabel <- nextLabel "while.cond"
  bodyLabel <- nextLabel "while.body"
  endLabel <- nextLabel "while.end"

  appendBrIfNotTerminated condLabel
  emitCurrentBlockAndStartNew condLabel
  codegenBoolExpr condExpr bodyLabel endLabel
  sealBlock bodyLabel
  sealBlock endLabel

  emitCurrentBlockAndStartNew bodyLabel
  codegenBlock block
  appendBrIfNotTerminated condLabel

  sealBlock condLabel
  emitCurrentBlockAndStartNew endLabel
codegenStmt (AST.SExpr e) = do
  codegenExpr e
  pure ()
codegenStmt (AST.SIncr (AST.Expr AST.ASTInt (AST.EVariable varName))) = do
  block <- gets cgCurrentBlock
  var <- readVariable L.Int varName block
  var' <- nextRegister IntR
  emit $ L.Add var' L.Int var "1"
  writeVariable varName block var'
codegenStmt (AST.SIncr lval@(AST.Expr AST.ASTInt _)) = do
  lvalAddr' <- codegenLVal lval
  lval' <- nextRegister IntR
  emit $ L.Load lval' L.Int lvalAddr'
  var' <- nextRegister IntR
  emit $ L.Add var' L.Int lval' "1"
  emit $ L.Store L.Int var' (L.Ptr L.Int) lvalAddr'
codegenStmt (AST.SDecr (AST.Expr AST.ASTInt (AST.EVariable varName))) = do
  block <- gets cgCurrentBlock
  var <- readVariable L.Int varName block
  var' <- nextRegister IntR
  emit $ L.Sub var' L.Int var "1"
  writeVariable varName block var'
codegenStmt (AST.SDecr lval@(AST.Expr AST.ASTInt _)) = do
  lvalAddr' <- codegenLVal lval
  lval' <- nextRegister IntR
  emit $ L.Load lval' L.Int lvalAddr'
  var' <- nextRegister IntR
  emit $ L.Sub var' L.Int lval' "1"
  emit $ L.Store L.Int var' (L.Ptr L.Int) lvalAddr'
codegenStmt (AST.SFor forVarT forVarName forCollection forBlock) = do
  let arrForVarT = AST.ASTArr forVarT
  idxName <- nextForIdx "i"
  idxSize <- nextForIdx "size"
  idxCollection <- nextForIdx "collection"
  let arrVar = AST.Expr arrForVarT $ AST.EVariable idxCollection
  let forStmt =
        AST.SBlock $
          AST.Block
            [ -- int i = 0;
              AST.SVarDecl idxName AST.ASTInt $ AST.Expr AST.ASTInt $ AST.EInt 0,
              -- T[] collection = forCollection;
              AST.SVarDecl idxCollection arrForVarT forCollection,
              -- int size = forCollection.size;
              AST.SVarDecl idxSize AST.ASTInt $ AST.Expr AST.ASTInt $ AST.EArrLength forVarT arrVar,
              -- while (i < size) {
              AST.SWhile
                ( AST.Expr
                    AST.ASTBool
                    ( AST.ELt
                        (AST.Expr AST.ASTInt (AST.EVariable idxName))
                        (AST.Expr AST.ASTInt (AST.EVariable idxSize))
                    )
                )
                $ AST.Block
                  [ -- forVarT forVarName = forCollection[i];
                    AST.SVarDecl forVarName forVarT $ AST.Expr forVarT $ AST.EArrGet arrVar $ AST.Expr AST.ASTInt $ AST.EVariable idxName,
                    -- forBlock;
                    AST.SBlock forBlock,
                    -- i = i + 1;
                    AST.SIncr $ AST.Expr AST.ASTInt $ AST.EVariable idxName
                  ]
            ]
  codegenStmt forStmt
codegenStmt stmt = do
  trace ("codegenStmt: " ++ show stmt) $ pure ()
  error "Not implemented"

-- book --

writeVariablesToBlock :: BlockID -> [(String, String)] -> CGM ()
writeVariablesToBlock blockId variables = do
  mapM_ go variables
  where
    go :: (String, String) -> CGM ()
    go (name, value) = do
      writeVariable name blockId value

getBlockPreds :: BlockID -> CGM [BlockID]
getBlockPreds blockId = do
  trace ("getBlockPreds: " ++ blockId) $ pure ()
  blockPreds <- gets cgBlockPreds
  trace ("getBlockPreds: " ++ blockId ++ " blockPreds:" ++ show blockPreds) $ pure ()
  pure $ Map.findWithDefault [] blockId blockPreds

getBlockPred :: Int -> BlockID -> CGM BlockID
getBlockPred i blockId = do
  blockPreds <- getBlockPreds blockId
  pure $ blockPreds !! i

newPhi :: L.Type -> BlockID -> CGM PhiID
newPhi ty blockId = newPhiWithOperands ty blockId []

newPhiWithOperands :: L.Type -> BlockID -> [(BlockID, String)] -> CGM PhiID
newPhiWithOperands ty blockId operands = do
  nextPhi <- nextPhi
  let phiObj = PhiObj nextPhi blockId operands ty
  modify $ \state -> state {cgPhis = Map.insert nextPhi phiObj $ cgPhis state}
  addPhiToBlock blockId nextPhi
  pure nextPhi

addPhiToBlock :: BlockID -> PhiID -> CGM ()
addPhiToBlock blockId phiId = do
  phis <- gets cgBlocksPhis
  let phis' = Map.insert blockId (phis Map.! blockId ++ [phiId]) phis
  modify $ \state -> state {cgBlocksPhis = phis'}

newUndef :: CGM String
newUndef = nextPhiPrefixed "undef.phi"

saveIncompletePhi :: String -> BlockID -> PhiID -> CGM ()
saveIncompletePhi variable blockId phiId = do
  trace ("saveIncompletePhi: " ++ variable ++ " " ++ blockId ++ " " ++ phiId) $ pure ()
  incompletePhis <- gets cgIncompletePhi
  let incompletePhi = Map.lookup blockId incompletePhis
  let newIncompletePhi'Block = case incompletePhi of
        Just incompletePhi'Block -> Map.insert variable phiId incompletePhi'Block
        Nothing -> Map.singleton variable phiId
  let newIncompletePhis = Map.insert blockId newIncompletePhi'Block incompletePhis
  modify $ \state -> state {cgIncompletePhi = newIncompletePhis}

writeVariable :: String -> BlockID -> String -> CGM ()
writeVariable variable block value = do
  trace ("writeVariable: " ++ variable ++ " " ++ block ++ " " ++ value) $ pure ()
  currentDef <- gets cgCurrentDef
  let blockMap = Map.findWithDefault Map.empty variable currentDef
  let newBlockMap = Map.insert block value blockMap
  let newCurrentDef = Map.insert variable newBlockMap currentDef
  modify $ \state -> state {cgCurrentDef = newCurrentDef}

readVariable :: L.Type -> String -> BlockID -> CGM String
readVariable ty variable block = do
  trace ("readVariable: " ++ variable ++ " " ++ block) $ pure ()
  currentDef <- gets cgCurrentDef
  let blockMap = Map.findWithDefault Map.empty variable currentDef
  let value = Map.lookup block blockMap
  case value of
    Just value' -> pure value'
    Nothing -> readVariableRecursive ty variable block

readVariableRecursive :: L.Type -> String -> BlockID -> CGM String
readVariableRecursive ty varName blockId = do
  trace ("readVariableRecursive: " ++ varName ++ " " ++ blockId) $ pure ()
  sealedBlocks <- gets cgSealedBlocks
  val <-
    if not $ Set.member blockId sealedBlocks
      then do
        trace ("readVariableRecursive: " ++ varName ++ " " ++ blockId ++ " not sealed") $ pure ()
        val <- newPhi ty blockId
        saveIncompletePhi varName blockId val
        pure val
      else do
        preds <- getBlockPreds blockId
        let predsSize = length preds
        if predsSize == 1
          then do
            trace ("readVariableRecursive: " ++ varName ++ " " ++ blockId ++ " predsSize == 1") $ pure ()
            readVariable ty varName $ head preds
          else do
            trace ("readVariableRecursive: " ++ varName ++ " " ++ blockId ++ " predsSize /= 1") $ pure ()
            val <- newPhi ty blockId
            writeVariable varName blockId val
            addPhiOperands varName val
  writeVariable varName blockId val
  pure val

sealBlock :: BlockID -> CGM ()
sealBlock blockId = do
  incompletePhis <- gets cgIncompletePhi
  let incompletePhi = Map.lookup blockId incompletePhis
  case incompletePhi of
    Just incompletePhi'Block -> do
      modify $ \state -> state {cgIncompletePhi = Map.delete blockId incompletePhis}
      mapM_ go $ Map.toList incompletePhi'Block
    Nothing -> pure ()
  modify $ \state -> state {cgSealedBlocks = Set.insert blockId $ cgSealedBlocks state}
  where
    go :: (String, PhiID) -> CGM ()
    go (variable, phiId) = do
      void $ addPhiOperands variable phiId

addPhiOperands :: String -> PhiID -> CGM String
addPhiOperands varName phiId = do
  trace ("addPhiOperands: " ++ varName ++ " " ++ phiId) $ pure ()
  phiObj <- getPhiObj phiId
  let blockId = phiBlock phiObj
  let ty = phiHeldType phiObj
  preds <- getBlockPreds blockId
  trace ("addPhiOperands: " ++ varName ++ " " ++ phiId ++ " blockID:" ++ show blockId) $ pure ()
  trace ("addPhiOperands: " ++ varName ++ " " ++ phiId ++ " preds:" ++ show preds) $ pure ()
  mapM_ (forEachPred ty phiId blockId) preds
  tryRemoveTrivialPhi phiId
  where
    forEachPred :: L.Type -> PhiID -> BlockID -> BlockID -> CGM ()
    forEachPred ty phiId blockId pred = do
      val <- readVariable ty varName pred
      phiObj <- getPhiObj phiId
      let phiObj' = PhiObj phiId blockId (phiOperands phiObj ++ [(pred, val)]) ty
      modify $ \state -> state {cgPhis = Map.insert phiId phiObj' $ cgPhis state}

data TrivialPhiCases
  = TPNone
  | TPReturn PhiID
  | TPOperand String

unwrapTrivialPhiCases :: TrivialPhiCases -> CGM String
unwrapTrivialPhiCases TPNone = newUndef
unwrapTrivialPhiCases (TPReturn phiId) = pure phiId
unwrapTrivialPhiCases (TPOperand op) = pure op

tryRemoveTrivialPhi :: PhiID -> CGM String
tryRemoveTrivialPhi phiId = do
  trace ("tryRemoveTrivialPhi: " ++ phiId) $ pure ()
  phiObj <- getPhiObj phiId
  let operands = phiOperands phiObj
  trace ("tryRemoveTrivialPhi: " ++ phiId ++ " operands: " ++ show operands) $ pure ()
  some <- unwrapTrivialPhiCases =<< foldM goPhiOperands TPNone operands
  if some == phiId
    then pure phiId
    else do
      phiOtherOperands <- allPhisThatUseMe phiId
      phiReplaceBy phiId some
      mapM_ tryRemoveTrivialPhi phiOtherOperands
      trace ("tryRemoveTrivialPhi: " ++ phiId ++ " some: " ++ some) $ pure ()
      pure some
  where
    goPhiOperands :: TrivialPhiCases -> (BlockID, String) -> CGM TrivialPhiCases
    goPhiOperands acc@(TPOperand some) (_, op) | some == op = pure acc
    goPhiOperands acc (_, op) | op == phiId = pure acc
    goPhiOperands TPOperand {} _ = pure $ TPReturn phiId
    goPhiOperands acc@TPReturn {} _ = pure acc
    goPhiOperands _ (_, op) = pure $ TPOperand op

allPhisThatUseMe :: PhiID -> CGM [PhiID]
allPhisThatUseMe myId = do
  phis <- gets cgPhis
  let phis' = Map.toList phis
  let phis'' = filter (\(_, phiObj) -> myId `elem` map snd (phiOperands phiObj)) phis'
  let phis''' = map fst phis''
  pure $ filter (/= myId) phis'''

-- # Reroute all uses of phi to same and remove phi
phiReplaceBy :: PhiID -> String -> CGM ()
phiReplaceBy phiId some = do
  runStatePhiReplacer phiId some
  modify $ \state@CGState {cgPhis = phis} -> state {cgPhis = Map.delete phiId phis}
  modify $ \state@CGState {cgBlocksPhis = blocksPhis} -> state {cgBlocksPhis = Map.map (filter (/= phiId)) blocksPhis}
  modify $ \state@CGState {cgIncompletePhi = incompletePhi} -> state {cgIncompletePhi = Map.map (Map.mapMaybeWithKey (\_ x -> if x == phiId then Nothing else Just x)) incompletePhi}

-- book --

sizeofType :: AST.ASTType -> CGM String
sizeofType AST.ASTInt = pure "4" -- i32
sizeofType AST.ASTBool = pure "1" -- i1
sizeofType AST.ASTStr = do
  sizePtrR <- nextRegister PtrR
  intSizeR <- nextRegister IntR
  let ty = L.UnsafeRaw "%struct.String"
  emit $ L.Gep sizePtrR ty "null" ["1"]
  emit $ L.Ptrtoint intSizeR L.String sizePtrR L.Int
  pure intSizeR
sizeofType (AST.ASTArr t) = do
  let arrT = L.Tuple [L.Int, L.Ptr $ L.fromAstType t]
  sizePtrR <- nextRegister PtrR
  intSizeR <- nextRegister IntR
  emit $ L.Gep sizePtrR arrT "null" ["1"]
  emit $ L.Ptrtoint intSizeR (L.Ptr arrT) sizePtrR L.Int
  pure intSizeR
sizeofType (AST.ASTClass name) = do
  let arrT = L.UnsafeRaw $ "%" ++ name
  sizePtrR <- nextRegister PtrR
  intSizeR <- nextRegister IntR
  emit $ L.Gep sizePtrR arrT "null" ["1"]
  emit $ L.Ptrtoint intSizeR (L.Ptr arrT) sizePtrR L.Int
  pure intSizeR
sizeofType AST.ASTVoid = error "sizeofType: void"

literalArrayToi8 :: Int -> String
literalArrayToi8 i = "[" ++ show i ++ " x i8]*"

codegenExpr :: AST.Expr -> CGM String
codegenExpr (AST.Expr ty (AST.EVariable name)) = do
  block <- gets cgCurrentBlock
  let ty' = L.fromAstType ty
  readVariable ty' name block
codegenExpr (AST.Expr AST.ASTInt (AST.EInt i)) =
  pure $ show i
codegenExpr (AST.Expr AST.ASTStr (AST.EString s)) = do
  let stringSize = length s
  strReg <- nextRegister StringR
  tempReg <- nextRegister PtrR
  strLiteralName <- addToStringPool s
  emit $ L.Bitcast tempReg (literalArrayToi8 stringSize) strLiteralName "i8*"
  emit $ L.Call strReg L.String "_new_str_from_literal" [(L.voidPtrType, tempReg), (L.Int, show stringSize)]
  pure strReg
codegenExpr (AST.Expr AST.ASTBool (AST.EBool b)) =
  pure $ if b then "true" else "false"
codegenExpr (AST.Expr AST.ASTInt (AST.EAdd e1 e2)) = do
  e1' <- codegenExpr e1
  e2' <- codegenExpr e2
  maybeRegister <- readLCSE LCSEAdd e1' e2'
  case maybeRegister of
    Just r -> pure r
    Nothing -> do
      r <- nextRegister IntR
      emit $ L.Add r L.Int e1' e2'
      saveLCSE LCSEAdd e1' e2' r
      pure r
codegenExpr (AST.Expr AST.ASTStr (AST.EAdd e1 e2)) = do
  e1' <- codegenExpr e1
  e2' <- codegenExpr e2
  r <- nextRegister StringR
  emit $ L.Call r L.String "_concat_str" [(L.String, e1'), (L.String, e2')]
  pure r
codegenExpr (AST.Expr AST.ASTInt (AST.ESub e1 e2)) = do
  e1' <- codegenExpr e1
  e2' <- codegenExpr e2
  maybeRegister <- readLCSE LCSESub e1' e2'
  case maybeRegister of
    Just r -> pure r
    Nothing -> do
      r <- nextRegister IntR
      emit $ L.Sub r L.Int e1' e2'
      saveLCSE LCSESub e1' e2' r
      pure r
codegenExpr (AST.Expr AST.ASTInt (AST.EMul e1 e2)) = do
  e1' <- codegenExpr e1
  e2' <- codegenExpr e2
  maybeRegister <- readLCSE LCSEMul e1' e2'
  case maybeRegister of
    Just r -> pure r
    Nothing -> do
      r <- nextRegister IntR
      emit $ L.Mul r L.Int e1' e2'
      saveLCSE LCSEMul e1' e2' r
      pure r
codegenExpr (AST.Expr AST.ASTInt (AST.EDiv e1 e2)) = do
  e1' <- codegenExpr e1
  e2' <- codegenExpr e2
  maybeRegister <- readLCSE LCSEDiv e1' e2'
  case maybeRegister of
    Just r -> pure r
    Nothing -> do
      r <- nextRegister IntR
      emit $ L.Div r L.Int e1' e2'
      saveLCSE LCSEDiv e1' e2' r
      pure r
codegenExpr (AST.Expr AST.ASTInt (AST.EMod e1 e2)) = do
  e1' <- codegenExpr e1
  e2' <- codegenExpr e2
  maybeRegister <- readLCSE LCSEMod e1' e2'
  case maybeRegister of
    Just r -> pure r
    Nothing -> do
      r <- nextRegister IntR
      emit $ L.Mod r L.Int e1' e2'
      saveLCSE LCSEMod e1' e2' r
      pure r
codegenExpr (AST.Expr AST.ASTInt (AST.ENeg e)) = do
  e' <- codegenExpr e
  r <- nextRegister IntR
  emit $ L.Neg r L.Int e'
  pure r
codegenExpr (AST.Expr ty (AST.EApplication fName args)) = do
  let ty' = L.fromAstType ty
  let tys = map (L.fromAstType . AST.exprType) args
  callTys <- getFunctionRealCallArgument fName
  args' <- mapM codegenExpr args
  args'' <- mapM castArgumentForCall (zip (zip tys args') callTys)
  let args''' = map snd args''
  r <- nextRegister $ registerTypeOfType ty'
  emit $ L.Call r ty' fName (zip callTys args''')
  pure r
codegenExpr (AST.Expr AST.ASTBool (AST.EEq e1@AST.Expr {AST.exprType = AST.ASTStr} e2)) = codegenStr e1 e2 "_str_eq_i1"
codegenExpr (AST.Expr AST.ASTBool (AST.ENe e1@AST.Expr {AST.exprType = AST.ASTStr} e2)) = codegenStr e1 e2 "_str_neq_i1"
codegenExpr (AST.Expr AST.ASTBool (AST.EEq e1 e2)) = do codegenIcmp e1 e2 L.Eq
codegenExpr (AST.Expr AST.ASTBool (AST.ENe e1 e2)) = do codegenIcmp e1 e2 L.Ne
codegenExpr (AST.Expr AST.ASTBool (AST.ELt e1 e2)) = do codegenIcmp e1 e2 L.Slt
codegenExpr (AST.Expr AST.ASTBool (AST.ELe e1 e2)) = do codegenIcmp e1 e2 L.Sle
codegenExpr (AST.Expr AST.ASTBool (AST.EGt e1 e2)) = do codegenIcmp e1 e2 L.Sgt
codegenExpr (AST.Expr AST.ASTBool (AST.EGe e1 e2)) = do codegenIcmp e1 e2 L.Sge
codegenExpr e@(AST.Expr AST.ASTBool (AST.EAnd e1 e2)) = do
  trueLabel <- nextLabel "and.true"
  falseLabel <- nextLabel "and.false"
  mergeLabel <- nextLabel "and.merge"

  codegenBoolExpr e trueLabel falseLabel
  sealBlock trueLabel
  sealBlock falseLabel

  emitCurrentBlockAndStartNew trueLabel
  appendBrIfNotTerminated mergeLabel
  emitCurrentBlockAndStartNew falseLabel
  appendBrIfNotTerminated mergeLabel

  sealBlock mergeLabel
  emitCurrentBlockAndStartNew mergeLabel
  newPhiWithOperands L.Bool mergeLabel [(trueLabel, "true"), (falseLabel, "false")]
codegenExpr e@(AST.Expr AST.ASTBool (AST.EOr e1 e2)) = do
  trueLabel <- nextLabel "or.true"
  falseLabel <- nextLabel "or.false"
  mergeLabel <- nextLabel "or.merge"

  codegenBoolExpr e trueLabel falseLabel
  sealBlock trueLabel
  sealBlock falseLabel

  emitCurrentBlockAndStartNew trueLabel
  appendBrIfNotTerminated mergeLabel
  emitCurrentBlockAndStartNew falseLabel
  appendBrIfNotTerminated mergeLabel

  sealBlock mergeLabel
  emitCurrentBlockAndStartNew mergeLabel
  newPhiWithOperands L.Bool mergeLabel [(trueLabel, "true"), (falseLabel, "false")]
codegenExpr (AST.Expr AST.ASTBool (AST.ENot e)) = do
  e' <- codegenExpr e
  r <- nextRegister BoolR
  emit $ L.Not r L.Bool e'
  pure r
codegenExpr (AST.Expr arrT@(AST.ASTArr arrDataT) (AST.ENewArray lenE)) = do
  singleElementSize <- sizeofType arrDataT
  arrSize <- sizeofType arrT
  let resultArrayType = L.fromAstType arrT
  let arrayDataType = L.Ptr $ L.fromAstType arrDataT
  let resultArrayTypeWithoutPtr = L.Tuple [L.Int, arrayDataType]
  lenE' <- codegenExpr lenE
  rawDataR <- nextRegister PtrR
  dataR <- nextRegister PtrR
  arrStructRawR <- nextRegister PtrR
  arrStructR <- nextRegister PtrR
  arrStructLenR <- nextRegister PtrR
  arrStructDataR <- nextRegister PtrR
  case arrDataT of
    AST.ASTStr -> do
      tempReg <- nextRegister PtrR
      emit $ L.Bitcast tempReg (literalArrayToi8 0) "@str.empty" "i8*"
      emit $ L.Call rawDataR L.voidPtrType "_memory_allocate_many_strings" [(L.Int, lenE'), (L.voidPtrType, tempReg)]
    _ -> emit $ L.Call rawDataR L.voidPtrType "_memory_allocate_many" [(L.Int, lenE'), (L.Int, singleElementSize)]
  emit $ L.BitcastT dataR L.voidPtrType rawDataR arrayDataType
  emit $ L.Call arrStructRawR L.voidPtrType "_memory_allocate" [(L.Int, arrSize)]
  emit $ L.BitcastT arrStructR L.voidPtrType arrStructRawR resultArrayType
  emit $ L.Gep arrStructLenR resultArrayTypeWithoutPtr arrStructR ["0", "0"]
  emit $ L.Store L.Int lenE' (L.Ptr L.Int) arrStructLenR
  emit $ L.Gep arrStructDataR resultArrayTypeWithoutPtr arrStructR ["0", "1"]
  emit $ L.Store arrayDataType dataR (L.Ptr arrayDataType) arrStructDataR
  pure arrStructR
codegenExpr e@(AST.Expr elemT (AST.EArrGet arrE elemE)) = do
  elemPtrR <- codegenLVal e
  r <- nextRegister PtrR
  emit $ L.Load r (L.fromAstType elemT) elemPtrR
  pure r
codegenExpr (AST.Expr AST.ASTInt (AST.EArrLength ty arrE)) = do
  let arrT = L.Tuple [L.Int, L.Ptr $ L.fromAstType ty]
  arrPtrR <- codegenExpr arrE
  ptrR <- nextRegister PtrR
  r <- nextRegister IntR
  emit $ L.Gep ptrR arrT arrPtrR ["0", "0"]
  emit $ L.Load r L.Int ptrR
  pure r
codegenExpr (AST.Expr _ AST.ENull) = do
  pure "null"
codegenExpr (AST.Expr (AST.ASTClass className) (AST.ENewObject name)) = do
  unless (className == name) $ error $ "codegenExpr: " ++ show className ++ " " ++ show name
  classObjR <- nextRegister PtrR
  emit $ L.Call classObjR (L.Class name) (name ++ ".constructor") []
  pure classObjR
codegenExpr e@(AST.Expr ty (AST.EFieldGet objE fName)) = do
  fieldR <- codegenLVal e
  r <- nextRegister PtrR
  emit $ L.Load r (L.fromAstType ty) fieldR
  pure r
codegenExpr (AST.Expr ty (AST.EMethodCall objE@AST.Expr {AST.exprType = AST.ASTClass className} fName args)) = do
  let ty' = L.fromAstType ty
  let tys = map (L.fromAstType . AST.exprType) args
  args' <- mapM codegenExpr args
  objE' <- codegenExpr objE
  let vTablePtrT = L.UnsafeRaw $ "%" ++ className ++ ".vtable"

  methodIdx <- calcMethodIndex className fName
  vTablePtrR <- nextRegister PtrR
  emit $ L.Gep vTablePtrR (L.UnsafeRaw $ "%" ++ className) objE' ["0", "0"]
  vTableR <- nextRegister PtrR
  emit $ L.Load vTableR (L.Ptr vTablePtrT) vTablePtrR
  fName <- nextRegister PtrR
  emit $ L.Gep fName vTablePtrT vTableR ["0", show methodIdx]

  mType@(L.Fun mRetType mArgs@(a : as)) <- getMethodType methodIdx className
  fNameCallable <- nextRegister PtrR
  emit $ L.Load fNameCallable mType fName

  let realArgs = (L.Class className, objE') : zip tys args'
  callableArgs <- mapM castArgumentForCall (zip realArgs mArgs)

  r <- nextRegister $ registerTypeOfType ty'
  emit $ L.CallM r ty' fNameCallable callableArgs
  pure r
codegenExpr e@(AST.Expr _ _) = do
  trace ("codegenExpr: " ++ show e) $ pure ()
  undefined

castArgumentForCall :: ((L.Type, String), L.Type) -> CGM (L.Type, String)
castArgumentForCall ((argType, arg), expectedType) = do
  if L.codegen argType /= L.codegen expectedType
    then do
      arg' <- nextRegister $ registerTypeOfType expectedType
      emit $ L.BitcastT arg' argType arg expectedType
      pure (expectedType, arg')
    else pure (expectedType, arg)

codegenIcmp :: AST.Expr -> AST.Expr -> L.Cond -> CGM String
codegenIcmp e1@AST.Expr {AST.exprType = t1} e2@AST.Expr {AST.exprType = t2} cond = do
  e1' <- codegenExpr e1
  e2'' <- codegenExpr e2
  e2' <-
    if t1 == t2
      then pure e2''
      else do
        r2 <- nextRegister $ registerTypeOfType $ L.fromAstType t1
        emit $ L.BitcastT r2 (L.fromAstType t2) e2'' (L.fromAstType t1)
        pure r2
  r <- nextRegister BoolR
  let ty = L.fromAstType t1
  emit $ L.Icmp r cond ty e1' e2'
  pure r

codegenStr :: AST.Expr -> AST.Expr -> String -> CGM String
codegenStr e1 e2 fName = do
  e1' <- codegenExpr e1
  e2' <- codegenExpr e2
  r <- nextRegister BoolR
  emit $ L.Call r L.Bool fName [(L.String, e1'), (L.String, e2')]
  pure r

codegenBoolExpr :: AST.Expr -> String -> String -> CGM ()
codegenBoolExpr (AST.Expr AST.ASTBool (AST.EAnd e1 e2)) trueBlock falseBlock = do
  lazyLabel <- nextLabel "and.lazy"
  codegenBoolExpr e1 lazyLabel falseBlock
  sealBlock lazyLabel
  emitCurrentBlockAndStartNew lazyLabel
  codegenBoolExpr e2 trueBlock falseBlock
codegenBoolExpr (AST.Expr AST.ASTBool (AST.EOr e1 e2)) trueBlock falseBlock = do
  lazyLabel <- nextLabel "or.lazy"
  codegenBoolExpr e1 trueBlock lazyLabel
  sealBlock lazyLabel
  emitCurrentBlockAndStartNew lazyLabel
  codegenBoolExpr e2 trueBlock falseBlock
codegenBoolExpr (AST.Expr AST.ASTBool (AST.ENot e)) trueBlock falseBlock = do
  codegenBoolExpr e falseBlock trueBlock
codegenBoolExpr e@(AST.Expr AST.ASTBool _) trueBlock falseBlock = do
  e' <- codegenExpr e
  appendCbrIfNotTerminated e' trueBlock falseBlock
codegenBoolExpr e t f = error $ "codegenBoolExpr: " ++ show e ++ " " ++ t ++ " " ++ f

codegenLVal :: AST.Expr -> CGM String
codegenLVal (AST.Expr t (AST.EArrGet eArr eInx)) = do
  let dataT = L.fromAstType t
  let dataPtrT = L.Ptr dataT
  eArr' <- codegenExpr eArr
  eInx' <- codegenExpr eInx
  dataR <- nextRegister PtrR
  arrR <- nextRegister PtrR
  r <- nextRegister PtrR
  let t' = L.Tuple [L.Int, dataPtrT]
  emit $ L.Gep dataR t' eArr' ["0", "1"]
  emit $ L.Load arrR dataPtrT dataR
  emit $ L.Gep r dataT arrR [eInx']
  pure r
codegenLVal e@(AST.Expr _ (AST.EFieldGet eObj@AST.Expr {AST.exprType = AST.ASTClass className} fName)) = do
  eObjR <- codegenExpr eObj
  r <- nextRegister PtrR
  idx <- calcFieldIndex className fName
  emit $ L.Gep r (L.UnsafeRaw $ "%" ++ className) eObjR ["0", show idx]
  pure r
codegenLVal e = do
  trace ("codegenLVal: " ++ show e) $ pure ()
  undefined

calcFieldIndex :: String -> String -> CGM Int
calcFieldIndex className fieldName = do
  fields <- gets (AST.classFields . (Map.! className) . AST.astEnvClasses . cgAstEnv)
  pure $ 1 + findIndex fieldName (map AST.classFieldName fields)

calcMethodIndex :: String -> String -> CGM Int
calcMethodIndex className methodName = do
  methods <- gets (AST.classMethods . (Map.! className) . AST.astEnvClasses . cgAstEnv)
  pure $ findIndex methodName (map (AST.funName . AST.methodFun) methods)

getMethodType :: Int -> String -> CGM L.Type
getMethodType idx className = do
  methods <- gets (L.progClasses . cgProgram)
  let classObj = findClass className methods
  let allTypes = L.classVTableType classObj
  pure $ allTypes !! idx

findIndex :: String -> [String] -> Int
findIndex name names = go names 0
  where
    go :: [String] -> Int -> Int
    go [] _ = error $ "findIndex: " ++ name ++ " " ++ show names
    go (x : xs) i | x == name = i
    go (_ : xs) i = go xs (i + 1)

findClass :: String -> [L.Class] -> L.Class
findClass name classes = go classes
  where
    go :: [L.Class] -> L.Class
    go [] = error $ "findClass: " ++ name ++ " " ++ show classes
    go (x : xs) | L.className x == name = x
    go (_ : xs) = go xs

getFunctionRealCallArgument :: String -> CGM [L.Type]
getFunctionRealCallArgument name = do
  astTypes <- gets (AST.funArgs . (Map.! name) . AST.astEnvFunctions . cgAstEnv)
  pure $ map L.fromAstType astTypes

tryRemoveReallyTrivialBlocks :: CGM ()
tryRemoveReallyTrivialBlocks = do
  blocks <- gets cgGeneratedBlocks
  trivialBlocks <- mapMaybeM isReallyTrivialBlock blocks
  trace ("tryRemoveReallyTrivialBlocks: " ++ show trivialBlocks) $ pure ()
  unless (checkForCycles trivialBlocks) $ replaceBlocks trivialBlocks
  where
    isReallyTrivialBlock :: BlockObj -> CGM (Maybe (String, String))
    isReallyTrivialBlock BlockObj {blockId = block} | endsWith block ".entry" = pure Nothing
    isReallyTrivialBlock (BlockObj myName instructions (L.Br replaceMyWith)) | all isComment instructions = do
      phis <- collectAllPhisThatDependOnLabel myName
      trace ("isReallyTrivialBlock: " ++ myName ++ " " ++ show phis) $ pure ()
      if null phis
        then pure $ Just (myName, replaceMyWith)
        else pure Nothing
    isReallyTrivialBlock _ = pure Nothing

endsWith :: Eq a => [a] -> [a] -> Bool
endsWith x y = y `isSuffixOf` x

collectAllPhisThatDependOnLabel :: BlockID -> CGM [PhiID]
collectAllPhisThatDependOnLabel label = do
  phis <- gets cgPhis
  let phis' = Map.toList phis
  let phis'' = filter (\(_, phiObj) -> label `elem` map fst (phiOperands phiObj)) phis'
  let phis''' = map fst phis''
  pure phis'''

isComment :: L.Instructions -> Bool
isComment (L.Comment _) = True
isComment _ = False

data VisitedState = Visiting | Visited deriving (Eq, Ord, Show)

type Visited = Map.Map String VisitedState

checkForCycles :: [(String, String)] -> Bool
checkForCycles edges = do
  let edgedDfs = dfs edges
  let nodes = map fst edges
  let (bs, _) = runState (mapM edgedDfs nodes) Map.empty
  or bs

dfs :: [(String, String)] -> String -> State Visited Bool
dfs edges' node = do
  visited <- get
  case Map.lookup node visited of
    Just Visiting -> pure True
    Just Visited -> pure False
    Nothing -> do
      put $ Map.insert node Visiting visited
      let edges'' = filter (\(from, _) -> from == node) edges'
      let edges''' = map snd edges''
      results <- mapM (dfs edges') edges'''
      put $ Map.insert node Visited visited
      pure $ or results

replaceBlocks :: [(String, String)] -> CGM ()
replaceBlocks edges' = do
  trace ("replaceBlocks: " ++ show edges') $ pure ()
  let edges = reduceEdges edges'
  trace ("replaceBlocks: " ++ show edges) $ pure ()
  blocks <- gets cgGeneratedBlocks
  let blocks' = map (replaceBlock $ Map.fromList edges) blocks
  let blocks'' = filter (\(BlockObj myName _ _) -> myName `notElem` map fst edges) blocks'
  modify $ \state -> state {cgGeneratedBlocks = blocks''}

replaceBlock :: Map.Map String String -> BlockObj -> BlockObj
replaceBlock edges blockObj@(BlockObj myName instructions terminator) = do
  let findReplacement x = Map.findWithDefault x x edges
  let (comments, newTerminator) =
        case terminator of
          L.Br replaceMyWith -> ([], L.Br $ findReplacement replaceMyWith)
          L.Cbr cond trueBlock falseBlock -> do
            let trueReplacement = findReplacement trueBlock
            let falseReplacement = findReplacement falseBlock
            if trueReplacement == falseReplacement
              then ([L.Comment "Cond check was skipped due to label reduction"], L.Br trueReplacement)
              else ([], L.Cbr cond trueReplacement falseReplacement)
          _ -> ([], terminator)
  BlockObj myName (instructions ++ comments) newTerminator

reduceEdges :: [(String, String)] -> [(String, String)]
reduceEdges edges = do
  let (_, y) = runState (mapM (reduceEdges' . fst) edges) $ Map.fromList edges
  Map.toList y
  where
    reduceEdges' :: String -> State (Map.Map String String) ()
    reduceEdges' edge = do
      last <- findLast edge
      modify $ Map.insert edge last
    findLast :: String -> State (Map.Map String String) String
    findLast edge = do
      visited <- get
      case Map.lookup edge visited of
        Nothing -> pure edge
        Just last -> findLast last
