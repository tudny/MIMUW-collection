-- File generated by the BNF Converter (bnfc 2.9.4.1).

{-# LANGUAGE CPP #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE LambdaCase #-}
#if __GLASGOW_HASKELL__ <= 708
{-# LANGUAGE OverlappingInstances #-}
#endif

-- | Pretty-printer for Src.

module Src.Jabba.Print where

import Prelude
  ( ($), (.)
  , Bool(..), (==), (<)
  , Int, Integer, Double, (+), (-), (*)
  , String, (++)
  , ShowS, showChar, showString
  , all, elem, foldr, id, map, null, replicate, shows, span
  )
import Data.Char ( Char, isSpace )
import qualified Src.Jabba.Abs

-- | The top-level printing method.

printTree :: Print a => a -> String
printTree = render . prt 0

type Doc = [ShowS] -> [ShowS]

doc :: ShowS -> Doc
doc = (:)

render :: Doc -> String
render d = rend 0 False (map ($ "") $ d []) ""
  where
  rend
    :: Int        -- ^ Indentation level.
    -> Bool       -- ^ Pending indentation to be output before next character?
    -> [String]
    -> ShowS
  rend i p = \case
      "["      :ts -> char '[' . rend i False ts
      "("      :ts -> char '(' . rend i False ts
      "{"      :ts -> onNewLine i     p . showChar   '{'  . new (i+1) ts
      "}" : ";":ts -> onNewLine (i-1) p . showString "};" . new (i-1) ts
      "}"      :ts -> onNewLine (i-1) p . showChar   '}'  . new (i-1) ts
      [";"]        -> char ';'
      ";"      :ts -> char ';' . new i ts
      t  : ts@(s:_) | closingOrPunctuation s
                   -> pending . showString t . rend i False ts
      t        :ts -> pending . space t      . rend i False ts
      []           -> id
    where
    -- Output character after pending indentation.
    char :: Char -> ShowS
    char c = pending . showChar c

    -- Output pending indentation.
    pending :: ShowS
    pending = if p then indent i else id

  -- Indentation (spaces) for given indentation level.
  indent :: Int -> ShowS
  indent i = replicateS (2*i) (showChar ' ')

  -- Continue rendering in new line with new indentation.
  new :: Int -> [String] -> ShowS
  new j ts = showChar '\n' . rend j True ts

  -- Make sure we are on a fresh line.
  onNewLine :: Int -> Bool -> ShowS
  onNewLine i p = (if p then id else showChar '\n') . indent i

  -- Separate given string from following text by a space (if needed).
  space :: String -> ShowS
  space t s =
    case (all isSpace t', null spc, null rest) of
      (True , _   , True ) -> []              -- remove trailing space
      (False, _   , True ) -> t'              -- remove trailing space
      (False, True, False) -> t' ++ ' ' : s   -- add space if none
      _                    -> t' ++ s
    where
      t'          = showString t []
      (spc, rest) = span isSpace s

  closingOrPunctuation :: String -> Bool
  closingOrPunctuation [c] = c `elem` closerOrPunct
  closingOrPunctuation _   = False

  closerOrPunct :: String
  closerOrPunct = ")],;"

parenth :: Doc -> Doc
parenth ss = doc (showChar '(') . ss . doc (showChar ')')

concatS :: [ShowS] -> ShowS
concatS = foldr (.) id

concatD :: [Doc] -> Doc
concatD = foldr (.) id

replicateS :: Int -> ShowS -> ShowS
replicateS n f = concatS (replicate n f)

-- | The printer class does the job.

class Print a where
  prt :: Int -> a -> Doc

instance {-# OVERLAPPABLE #-} Print a => Print [a] where
  prt i = concatD . map (prt i)

instance Print Char where
  prt _ c = doc (showChar '\'' . mkEsc '\'' c . showChar '\'')

instance Print String where
  prt _ = printString

printString :: String -> Doc
printString s = doc (showChar '"' . concatS (map (mkEsc '"') s) . showChar '"')

mkEsc :: Char -> Char -> ShowS
mkEsc q = \case
  s | s == q -> showChar '\\' . showChar s
  '\\' -> showString "\\\\"
  '\n' -> showString "\\n"
  '\t' -> showString "\\t"
  s -> showChar s

prPrec :: Int -> Int -> Doc -> Doc
prPrec i j = if j < i then parenth else id

instance Print Integer where
  prt _ x = doc (shows x)

instance Print Double where
  prt _ x = doc (shows x)

instance Print Src.Jabba.Abs.Ident where
  prt _ (Src.Jabba.Abs.Ident i) = doc $ showString i
instance Print (Src.Jabba.Abs.Program' a) where
  prt i = \case
    Src.Jabba.Abs.PProgram _ instrs -> prPrec i 0 (concatD [prt 0 instrs])

instance Print (Src.Jabba.Abs.Instr' a) where
  prt i = \case
    Src.Jabba.Abs.DFun _ id_ args type_ block -> prPrec i 0 (concatD [doc (showString "fun"), prt 0 id_, doc (showString "("), prt 0 args, doc (showString ")"), doc (showString ":"), prt 0 type_, prt 0 block])
    Src.Jabba.Abs.DFunUnit _ id_ args block -> prPrec i 0 (concatD [doc (showString "fun"), prt 0 id_, doc (showString "("), prt 0 args, doc (showString ")"), prt 0 block])
    Src.Jabba.Abs.IUnit _ -> prPrec i 0 (concatD [doc (showString ";")])
    Src.Jabba.Abs.IIncr _ id_ -> prPrec i 0 (concatD [prt 0 id_, doc (showString "++"), doc (showString ";")])
    Src.Jabba.Abs.IDecr _ id_ -> prPrec i 0 (concatD [prt 0 id_, doc (showString "--"), doc (showString ";")])
    Src.Jabba.Abs.IAss _ id_ expr -> prPrec i 0 (concatD [prt 0 id_, doc (showString "="), prt 0 expr, doc (showString ";")])
    Src.Jabba.Abs.IRet _ expr -> prPrec i 0 (concatD [doc (showString "return"), prt 0 expr, doc (showString ";")])
    Src.Jabba.Abs.IRetUnit _ -> prPrec i 0 (concatD [doc (showString "return"), doc (showString ";")])
    Src.Jabba.Abs.IBreak _ -> prPrec i 0 (concatD [doc (showString "break"), doc (showString ";")])
    Src.Jabba.Abs.ICont _ -> prPrec i 0 (concatD [doc (showString "continue"), doc (showString ";")])
    Src.Jabba.Abs.IIf _ ifstmt -> prPrec i 0 (concatD [prt 0 ifstmt])
    Src.Jabba.Abs.IWhile _ expr block -> prPrec i 0 (concatD [doc (showString "while"), prt 0 expr, prt 0 block])
    Src.Jabba.Abs.IWhileFin _ expr block1 block2 -> prPrec i 0 (concatD [doc (showString "while"), prt 0 expr, prt 0 block1, doc (showString "finally"), prt 0 block2])
    Src.Jabba.Abs.IFor _ id_ expr1 expr2 block -> prPrec i 0 (concatD [doc (showString "for"), prt 0 id_, doc (showString "="), prt 0 expr1, doc (showString ".."), prt 0 expr2, prt 0 block])
    Src.Jabba.Abs.IExpr _ expr -> prPrec i 0 (concatD [prt 0 expr, doc (showString ";")])
    Src.Jabba.Abs.IDecl _ decl -> prPrec i 0 (concatD [prt 0 decl, doc (showString ";")])
    Src.Jabba.Abs.IBBlock _ block -> prPrec i 0 (concatD [prt 0 block])
    Src.Jabba.Abs.ITabAss _ id_ expr1 expr2 -> prPrec i 0 (concatD [prt 0 id_, doc (showString "["), prt 0 expr1, doc (showString "]"), doc (showString "="), prt 0 expr2, doc (showString ";")])

instance Print (Src.Jabba.Abs.Arg' a) where
  prt i = \case
    Src.Jabba.Abs.RefMutArg _ id_ type_ -> prPrec i 0 (concatD [doc (showString "var"), prt 0 id_, doc (showString ":"), prt 0 type_])
    Src.Jabba.Abs.RefConstArg _ id_ type_ -> prPrec i 0 (concatD [doc (showString "val"), prt 0 id_, doc (showString ":"), prt 0 type_])
    Src.Jabba.Abs.CopyMutArg _ id_ type_ -> prPrec i 0 (concatD [doc (showString "var"), prt 0 id_, doc (showString ":"), doc (showString "new"), prt 0 type_])
    Src.Jabba.Abs.CopyConstArg _ id_ type_ -> prPrec i 0 (concatD [doc (showString "val"), prt 0 id_, doc (showString ":"), doc (showString "new"), prt 0 type_])

instance Print [Src.Jabba.Abs.Arg' a] where
  prt _ [] = concatD []
  prt _ [x] = concatD [prt 0 x]
  prt _ (x:xs) = concatD [prt 0 x, doc (showString ","), prt 0 xs]

instance Print (Src.Jabba.Abs.Item' a) where
  prt i = \case
    Src.Jabba.Abs.DItemVal _ id_ type_ expr -> prPrec i 0 (concatD [prt 0 id_, doc (showString ":"), prt 0 type_, doc (showString "="), prt 0 expr])
    Src.Jabba.Abs.DItemAuto _ id_ expr -> prPrec i 0 (concatD [prt 0 id_, doc (showString "="), prt 0 expr])
    Src.Jabba.Abs.DItem _ id_ type_ -> prPrec i 0 (concatD [prt 0 id_, doc (showString ":"), prt 0 type_])

instance Print [Src.Jabba.Abs.Item' a] where
  prt _ [] = concatD []
  prt _ [x] = concatD [prt 0 x]
  prt _ (x:xs) = concatD [prt 0 x, doc (showString ","), prt 0 xs]

instance Print (Src.Jabba.Abs.Decl' a) where
  prt i = \case
    Src.Jabba.Abs.DVar _ items -> prPrec i 0 (concatD [doc (showString "var"), prt 0 items])
    Src.Jabba.Abs.DVal _ items -> prPrec i 0 (concatD [doc (showString "val"), prt 0 items])

instance Print (Src.Jabba.Abs.Block' a) where
  prt i = \case
    Src.Jabba.Abs.IBlock _ instrs -> prPrec i 0 (concatD [doc (showString "{"), prt 0 instrs, doc (showString "}")])

instance Print [Src.Jabba.Abs.Instr' a] where
  prt _ [] = concatD []
  prt _ (x:xs) = concatD [prt 0 x, prt 0 xs]

instance Print (Src.Jabba.Abs.IfStmt' a) where
  prt i = \case
    Src.Jabba.Abs.IfIf _ expr block -> prPrec i 0 (concatD [doc (showString "if"), prt 0 expr, prt 0 block])
    Src.Jabba.Abs.IfElse _ expr block1 block2 -> prPrec i 0 (concatD [doc (showString "if"), prt 0 expr, prt 0 block1, doc (showString "else"), prt 0 block2])
    Src.Jabba.Abs.IfElseIf _ expr block ifstmt -> prPrec i 0 (concatD [doc (showString "if"), prt 0 expr, prt 0 block, doc (showString "else"), prt 0 ifstmt])

instance Print (Src.Jabba.Abs.Type' a) where
  prt i = \case
    Src.Jabba.Abs.TInt _ -> prPrec i 0 (concatD [doc (showString "Integer")])
    Src.Jabba.Abs.TBool _ -> prPrec i 0 (concatD [doc (showString "Boolean")])
    Src.Jabba.Abs.TString _ -> prPrec i 0 (concatD [doc (showString "String")])
    Src.Jabba.Abs.TVoid _ -> prPrec i 0 (concatD [doc (showString "Unit")])
    Src.Jabba.Abs.TFun _ targs type_ -> prPrec i 0 (concatD [doc (showString "("), prt 0 targs, doc (showString ")"), doc (showString "->"), prt 0 type_])
    Src.Jabba.Abs.TTab _ type_ -> prPrec i 0 (concatD [doc (showString "["), prt 0 type_, doc (showString "]")])

instance Print (Src.Jabba.Abs.TArg' a) where
  prt i = \case
    Src.Jabba.Abs.TRefMutArg _ type_ -> prPrec i 0 (concatD [doc (showString "var"), prt 0 type_])
    Src.Jabba.Abs.TRefConstArg _ type_ -> prPrec i 0 (concatD [doc (showString "val"), prt 0 type_])
    Src.Jabba.Abs.TCopyMutArg _ type_ -> prPrec i 0 (concatD [doc (showString "var"), doc (showString "$"), prt 0 type_])
    Src.Jabba.Abs.TCopyConstArg _ type_ -> prPrec i 0 (concatD [doc (showString "val"), doc (showString "$"), prt 0 type_])

instance Print [Src.Jabba.Abs.TArg' a] where
  prt _ [] = concatD []
  prt _ [x] = concatD [prt 0 x]
  prt _ (x:xs) = concatD [prt 0 x, doc (showString ","), prt 0 xs]

instance Print (Src.Jabba.Abs.PlsOp' a) where
  prt i = \case
    Src.Jabba.Abs.OPlus _ -> prPrec i 0 (concatD [doc (showString "+")])
    Src.Jabba.Abs.OMinus _ -> prPrec i 0 (concatD [doc (showString "-")])

instance Print (Src.Jabba.Abs.MulOp' a) where
  prt i = \case
    Src.Jabba.Abs.OMul _ -> prPrec i 0 (concatD [doc (showString "*")])
    Src.Jabba.Abs.ODiv _ -> prPrec i 0 (concatD [doc (showString "/")])
    Src.Jabba.Abs.OMod _ -> prPrec i 0 (concatD [doc (showString "%")])

instance Print (Src.Jabba.Abs.NotOp' a) where
  prt i = \case
    Src.Jabba.Abs.ONot _ -> prPrec i 0 (concatD [doc (showString "!")])

instance Print (Src.Jabba.Abs.NegOp' a) where
  prt i = \case
    Src.Jabba.Abs.ONeg _ -> prPrec i 0 (concatD [doc (showString "-")])

instance Print (Src.Jabba.Abs.AndOp' a) where
  prt i = \case
    Src.Jabba.Abs.OAnd _ -> prPrec i 0 (concatD [doc (showString "&&")])

instance Print (Src.Jabba.Abs.OrOp' a) where
  prt i = \case
    Src.Jabba.Abs.OOr _ -> prPrec i 0 (concatD [doc (showString "||")])

instance Print (Src.Jabba.Abs.RelOp' a) where
  prt i = \case
    Src.Jabba.Abs.REq _ -> prPrec i 0 (concatD [doc (showString "==")])
    Src.Jabba.Abs.RNeq _ -> prPrec i 0 (concatD [doc (showString "!=")])
    Src.Jabba.Abs.RLt _ -> prPrec i 0 (concatD [doc (showString "<")])
    Src.Jabba.Abs.RGt _ -> prPrec i 0 (concatD [doc (showString ">")])
    Src.Jabba.Abs.RLeq _ -> prPrec i 0 (concatD [doc (showString "<=")])
    Src.Jabba.Abs.RGeq _ -> prPrec i 0 (concatD [doc (showString ">=")])

instance Print (Src.Jabba.Abs.Expr' a) where
  prt i = \case
    Src.Jabba.Abs.ITabAcc _ id_ expr -> prPrec i 6 (concatD [prt 0 id_, doc (showString "["), prt 0 expr, doc (showString "]")])
    Src.Jabba.Abs.ITabInit _ expr1 expr2 -> prPrec i 0 (concatD [doc (showString "["), prt 0 expr1, doc (showString ";"), prt 0 expr2, doc (showString "]")])
    Src.Jabba.Abs.ITabInitEls _ exprs -> prPrec i 0 (concatD [doc (showString "["), prt 0 exprs, doc (showString "]")])
    Src.Jabba.Abs.EVarName _ id_ -> prPrec i 6 (concatD [prt 0 id_])
    Src.Jabba.Abs.EIntLit _ n -> prPrec i 6 (concatD [prt 0 n])
    Src.Jabba.Abs.EBoolLitTrue _ -> prPrec i 6 (concatD [doc (showString "true")])
    Src.Jabba.Abs.EBoolLitFalse _ -> prPrec i 6 (concatD [doc (showString "false")])
    Src.Jabba.Abs.EUnitLiteral _ -> prPrec i 6 (concatD [doc (showString "unit")])
    Src.Jabba.Abs.EStringLit _ str -> prPrec i 6 (concatD [printString str])
    Src.Jabba.Abs.ERun _ expr exprs -> prPrec i 5 (concatD [prt 5 expr, doc (showString "("), prt 0 exprs, doc (showString ")")])
    Src.Jabba.Abs.ELambda _ args block -> prPrec i 5 (concatD [doc (showString "|"), prt 0 args, doc (showString "|"), prt 0 block])
    Src.Jabba.Abs.ELambdaEmpty _ block -> prPrec i 5 (concatD [doc (showString "||"), prt 0 block])
    Src.Jabba.Abs.ENeg _ negop expr -> prPrec i 5 (concatD [prt 0 negop, prt 6 expr])
    Src.Jabba.Abs.ENot _ notop expr -> prPrec i 5 (concatD [prt 0 notop, prt 6 expr])
    Src.Jabba.Abs.EMul _ expr1 mulop expr2 -> prPrec i 4 (concatD [prt 4 expr1, prt 0 mulop, prt 5 expr2])
    Src.Jabba.Abs.ESum _ expr1 plsop expr2 -> prPrec i 3 (concatD [prt 3 expr1, prt 0 plsop, prt 4 expr2])
    Src.Jabba.Abs.ERel _ expr1 relop expr2 -> prPrec i 2 (concatD [prt 2 expr1, prt 0 relop, prt 3 expr2])
    Src.Jabba.Abs.EBAnd _ expr1 andop expr2 -> prPrec i 1 (concatD [prt 2 expr1, prt 0 andop, prt 1 expr2])
    Src.Jabba.Abs.EBOr _ expr1 orop expr2 -> prPrec i 0 (concatD [prt 1 expr1, prt 0 orop, prt 0 expr2])
    Src.Jabba.Abs.ETer _ expr1 expr2 expr3 -> prPrec i 0 (concatD [prt 1 expr1, doc (showString "?"), prt 1 expr2, doc (showString ":"), prt 0 expr3])
    Src.Jabba.Abs.ELambdaExpr _ args expr -> prPrec i 0 (concatD [doc (showString "|"), prt 0 args, doc (showString "|"), prt 1 expr])
    Src.Jabba.Abs.ELambdaEmptEpr _ expr -> prPrec i 0 (concatD [doc (showString "||"), prt 1 expr])

instance Print [Src.Jabba.Abs.Expr' a] where
  prt _ [] = concatD []
  prt _ [x] = concatD [prt 0 x]
  prt _ (x:xs) = concatD [prt 0 x, doc (showString ","), prt 0 xs]
