import Control.Monad.Reader (ReaderT (runReaderT), MonadReader (ask, local), Reader, runReader)
import qualified Data.Map as Map
import Control.Monad.Except (ExceptT, MonadError (throwError), runExceptT)
import Control.Monad.State (StateT(runStateT))


type Var = String
data Exp = EInt Int
    | EOp  Op Exp Exp
    | EVar Var
    | ELet Var Exp Exp

data Op = OpAdd | OpMul | OpSub

type Env = Map.Map Var Int

resolveOp :: Op -> Int -> Int -> Int
resolveOp OpAdd = (+)
resolveOp OpMul = (*)
resolveOp OpSub = (-)


type IM a = ReaderT Env (Either String) a


evalExp :: Exp -> Either String Int
evalExp e = let result = runReaderT (go e) Map.empty in result
    where
        go :: Exp -> IM Int
        go (EInt n) = pure n
        go (EOp op e1 e2) = do
            e1' <- go e1
            e2' <- go e2
            let op' = resolveOp op
            pure $ op' e1' e2'
        go (EVar v) = do
            env <- ask
            case Map.lookup v env of 
                Nothing -> throwError "Value not defined"
                Just x -> pure x
        go (ELet v e1 e2) = do
            e1' <- go e1
            local (Map.insert v e1') $ go e2


sampleExp :: Exp
sampleExp = ELet "x" (EOp OpAdd (EInt 1) (EInt 2)) (EOp OpMul (EVar "x") (EInt 3))

sampleErrExp :: Exp
sampleErrExp = ELet "x" (EOp OpAdd (EInt 1) (EInt 2)) (EOp OpMul (EVar "y") (EInt 3))

-- >>> evalExp sampleExp
-- >>> evalExp sampleErrExp
-- Right 9
-- Left "Value not defined"



type IME a = ExceptT String (Reader Env) a

evalExp2 :: Exp -> Either String Int
evalExp2 e = let result = runReader (runExceptT (go e)) Map.empty in result
    where
        go :: Exp -> IME Int
        go (EInt n) = pure n
        go (EOp op e1 e2) = do
            n1 <- go e1
            n2 <- go e2
            pure $ resolveOp op n1 n2
        go (EVar v) = do
            env <- ask
            case Map.lookup v env of 
                Nothing -> throwError "Value not defined"
                Just x -> pure x
        go (ELet v e1 e2) = do
            e1' <- go e1
            local (Map.insert v e1') $ go e2

-- >>> evalExp2 sampleExp
-- >>> evalExp2 sampleErrExp
-- Right 9
-- Left "Value not defined"


type IMM a = ExceptT String ((->) Env) a

evalExp3 :: Exp -> Either String Int
evalExp3 e = let result = runExceptT (go e) Map.empty in result
    where
        go :: Exp -> IMM Int
        go (EInt n) = pure n
        go (EOp op e1 e2) = do
            n1 <- go e1
            n2 <- go e2
            pure $ resolveOp op n1 n2
        go (EVar v) = do
            env <- ask
            case Map.lookup v env of 
                Nothing -> throwError "Value not defined"
                Just x -> pure x
        go (ELet v e1 e2) = do
            e1' <- go e1
            local (Map.insert v e1') $ go e2

-- >>> evalExp3 sampleExp
-- >>> evalExp3 sampleErrExp
-- Right 9
-- Left "Value not defined"
