module AbstractSyntaxTree where

import Control.Monad.State
import Data.Map (Map, empty, insert, member)
import Instant.Abs

data Variables = Variables {getVars :: Map String Integer, getCount :: Integer}

data AST = AST {getStmts :: [IStmt], getVariables :: Variables}

type Err = Either String

type AM = StateT Variables Err

data IStmt
  = IAss String EExpr
  | IExpr EExpr

data EExpr'
  = EAdd EExpr EExpr
  | ESub EExpr EExpr
  | EMul EExpr EExpr
  | EDiv EExpr EExpr
  | EInt Integer
  | EVar String

data EExpr = EExpr' EExpr' Int Bool

parseTree :: Program -> Err AST
parseTree (Prog _ stmts) = do
  (stmts', vars) <- runStateT go (Variables empty 0)
  pure $ AST (reverse stmts') vars
  where
    go = foldM parseStmt' [] stmts

parseStmt' :: [IStmt] -> Stmt -> AM [IStmt]
parseStmt' stmts stmt = do
  stmt' <- parseStmt stmt
  pure (stmt' : stmts)

parseStmt :: Stmt -> AM IStmt
parseStmt (SAss _ (Ident ident) expr) = do
  e <- parseExpr expr
  vars <- get
  put $ computeIfAbsent vars ident
  pure $ IAss ident e
parseStmt (SExp _ expr) = do
  e <- parseExpr expr
  pure $ IExpr e

parseExpr :: Exp -> AM EExpr
parseExpr (ExpAdd _ e1 e2) = makeTwoExprExpr EAdd e1 e2
parseExpr (ExpSub _ e1 e2) = makeTwoExprExpr ESub e1 e2
parseExpr (ExpMul _ e1 e2) = makeTwoExprExpr EMul e1 e2
parseExpr (ExpDiv _ e1 e2) = makeTwoExprExpr EDiv e1 e2
parseExpr (ExpLit _ i) = pure $ EExpr' (EInt i) 1 False
parseExpr (ExpVar _ (Ident vname)) = do
  vars <- get
  if member vname (getVars vars)
    then pure $ EExpr' (EVar vname) 1 False
    else error $ "Variable " ++ vname ++ " not declared"

makeTwoExprExpr :: (EExpr -> EExpr -> EExpr') -> Exp -> Exp -> AM EExpr
makeTwoExprExpr f e1 e2 = do
  e1'@(EExpr' _ s1' _) <- parseExpr e1
  e2'@(EExpr' _ s2' _) <- parseExpr e2
  let doSwap' = s1' < s2'
  pure $
    EExpr'
      (f e1' e2')
      ( if doSwap'
          then max s2' (s1' + 1)
          else max s1' (s2' + 1)
      )
      doSwap'

computeIfAbsent :: Variables -> String -> Variables
computeIfAbsent vars name =
  if member name (getVars vars)
    then vars
    else
      let newCount = getCount vars + 1
       in Variables (insert name newCount (getVars vars)) newCount
