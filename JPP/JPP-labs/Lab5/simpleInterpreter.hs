import Control.Monad.State (State, runState, MonadState (get, put), gets, modify)
import qualified Data.Map as Map
import Control.Monad.Except


-- b. Korzystając z monady State, 
--  napisz interpreter języka Tiny 
--  (patrz przedmiot Semantyka i Weryfikacja Programów)

type Env = Map.Map Var Int

type Var = String
data Exp = EVar Var 
         | EInt Int 
         | EAdd Exp Exp 
         | ESub Exp Exp 
         | EMul Exp Exp 
         | EDiv Exp Exp
         | ELet Var Exp Exp
data BExp = BEq Exp Exp | BNot BExp
data Stmt = SSkip | SAss Var Exp | SComp Stmt Stmt
            | SIf BExp Stmt Stmt | SWhile BExp Stmt
            | SExp Exp

type IM a = ExceptT String (State Env) a

execStmt :: Env -> Stmt -> IO ()
execStmt e s = let (result, env) = runState (runExceptT (goS s)) e in do
    print result
    print env
    where
        goE :: Exp -> IM Int
        goB :: BExp -> IM Bool
        goS :: Stmt -> IM Int
        goS SSkip = do pure 0
        goS (SAss v e) = do
            n <- goE e
            modify $ Map.insert v n
            pure n  -- return n := E[|e|]
        goS (SComp s1 s2) = do
            s1' <- goS s1
            goS s2
        goS (SIf be s1 s2) = do
            b <- goB be
            goS (if b then s1 else s2)
        goS sw@(SWhile be s) = do
            b <- goB be
            if b then do
                goS s
                goS sw 
            else 
                pure 0
        goS (SExp e) = do goE e
        goB (BEq e1 e2) = do
            e1' <- goE e1
            e2' <- goE e2
            pure $ e1' == e2'
        goB (BNot be) = do
            be' <- goB be
            pure $ not be'
        goE (EVar v) = do 
            env <- get
            case Map.lookup v env of 
                Nothing -> throwError $ "Value " ++ v ++ " not defined"
                Just x -> pure x
        goE (EInt n) = do pure n
        goE (EAdd e1 e2) = do
            e1' <- goE e1
            e2' <- goE e2
            pure $ e1' + e2'
        goE (ESub e1 e2) = do
            e1' <- goE e1
            e2' <- goE e2
            pure $ e1' - e2'
        goE (EMul e1 e2) = do
            e1' <- goE e1
            e2' <- goE e2
            pure $ e1' * e2'
        goE (EDiv e1 e2) = do
            e1' <- goE e1
            e2' <- goE e2
            case e2' of
                0 -> throwError "Cannot divide by 0"
                n -> pure n
        goE (ELet v e1 e2) = do
            e1' <- goE e1
            saveEnv <- get
            modify $ Map.insert v e1'
            n <- goE e2
            put saveEnv
            pure n
    


sampleStmt :: Stmt
sampleStmt =
    SComp (
        SAss "fst" (EInt 0)
    ) (
        SComp (
            SAss "snd" (EInt 1)
        )
        (
            SComp (
                SAss "i" (EInt 0)
            )
            (
                SComp (
                    SWhile (
                        BNot (
                            BEq (EVar "i") (EVar "n")
                        )
                    )
                    (
                        SComp (
                            SAss "temp" (EVar "snd")
                        )
                        (
                            SComp (
                                SAss "snd" (EAdd (EVar "snd") (EVar "fst"))
                            )
                            (
                                SComp (
                                    SAss "fst" (EVar "temp")
                                )
                                (
                                    SAss "i" (EAdd (EVar "i") (EInt 1))
                                )
                            )
                        )
                    )
                )
                (
                    SExp (EVar "snd")
                )
            )
        )
    )


sampleProg :: Stmt
sampleProg = SComp (
        SAss "y" (
            ELet "x" (EInt 1) (EAdd (EVar "x") (EInt 1)) 
        )
    ) (
        SExp (EVar "x")
    )

-- y = let x = 1 in x + 1
-- x

-- >>> execStmt Map.empty sampleProg


main :: IO ()
main = do
    -- n <- readLn :: IO Int
    -- execStmt (Map.fromList [("n", n)]) sampleStmt
    execStmt Map.empty sampleProg

-- fst = 0
-- snd = 1
-- i = 0
-- while !(i == n) {
--     temp = snd
--     snd = snd + fst
--     fst = temp
--     i = i + 1
-- } 
-- snd  // <- return nth fibonacci number
