import Control.Monad.State (State, runState, MonadState (get, put), gets, modify)
import qualified Data.Map as Map
import Debug.Trace (trace)


-- b. Korzystając z monady State, 
--  napisz interpreter języka Tiny 
--  (patrz przedmiot Semantyka i Weryfikacja Programów)

type IntState = Map.Map Var Int

type Var = String
data Exp = EVar Var | EInt Int | EAdd Exp Exp | ESub Exp Exp
data BExp = BEq Exp Exp | BNot BExp
data Stmt = SSkip | SAss Var Exp | SComp Stmt Stmt
            | SIf BExp Stmt Stmt | SWhile BExp Stmt
            | SExp Exp

type IM a = State IntState a

execStmt :: IntState -> Stmt -> IO ()
execStmt e s = let (val, state) = runState (goS s) e in do
    print state
    print val
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
        goE (EVar v) = do gets $ Map.findWithDefault 0 v
        goE (EInt n) = do pure n
        goE (EAdd e1 e2) = do
            e1' <- goE e1
            e2' <- goE e2
            pure $ e1' + e2'
        goE (ESub e1 e2) = do
            e1' <- goE e1
            e2' <- goE e2
            pure $ e1' - e2'


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


main :: IO ()
main = do
    n <- readLn :: IO Int
    execStmt (Map.fromList [("n", n)]) sampleStmt

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
