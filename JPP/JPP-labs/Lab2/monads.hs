module Main where
import Data.Text.Internal.Fusion.Size (exactSize)
import Language.Haskell.TH (fromE)

-- 1. Rozważmy typ drzew trochę inny niż na wykladzie
data Tree a = Empty 
    | Node a (Tree a) (Tree a)

instance Show a => Show (Tree a) where
    show Empty = ""
    show (Node a l r) = "(" ++ show l ++ ") " ++ show a ++ " (" ++ show r ++ ")"


exampleTree = Node 1 (Node 2 Empty Empty) (Node 3 Empty Empty)
-- >>> show exampleTree
-- "(() 2 ()) 1 (() 3 ())"

instance Eq a => Eq (Tree a) where
    t1 == t2 = case (t1, t2) of
        (Empty, Empty) -> True
        (Node a1 l1 r1, Node a2 l2 r2) -> a1 == a2 && l1 == l2 && r1 == r2
        (_, _) -> False

-- fold from right
foldTree :: (a -> b -> b) -> b -> Tree a -> b
foldTree m v t = case t of 
    Empty -> v
    (Node a l r) -> foldTree m (m a (foldTree m v r)) l

toList :: Tree a -> [a]
toList = foldTree (:) []

-- >>> toList exampleTree
-- [2,1,3]


-- 2. Zdefniować data MyMaybe a = MyNothing | MyJust a oraz showsPrec dla niego tak aby nawiasy były tam gdzie trzeba, czyli np
data MyMaybe a = MyNothing | MyJust a

-- 3. Niech typ
-- newtype OrderedList a = OL [a]
-- reprezentuje listy uporządkowane niemalejąco

-- a. Uzupełnij instancje

-- instance  Ord a => Semigroup (OrderedList a) where
-- instance Ord a => Monoid (OrderedList a) where

newtype OrderedList a = OL [a]

instance  Ord a => Semigroup (OrderedList a) where
  (OL t1) <> (OL t2) = OL (merge t1 t2)
    where 
        merge :: Ord a => [a] -> [a] -> [a]
        merge [] t = t
        merge t [] = t
        merge (x1:xs1) (x2:xs2) 
            | x1 < x2   = x1 : merge xs1 (x2:xs2)
            | otherwise = x2 : merge (x1:xs1) xs2

l1 = OL [1, 2, 3, 6]
l2 = OL [4, 5, 7]

-- >>> let OL t = l1 <> l2 in t
-- [1,2,3,4,5,6,7]

-- b. Napisz funkcję eliminującą duplikaty (analogicznie do nub)
-- nubOrdered :: Ord a => OrderedList a -> OrderedList a

nubOrdered :: Ord a => OrderedList a -> OrderedList a
nubOrdered (OL t) = OL $ foldr (\ x a -> case a of { [] -> [x]; (s:ss) -> if x == s then a else x : a }) [] t 

-- >>> let OL t = nubOrdered (OL [0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 5, 6, 7, 7]) in t
-- [0,1,2,3,4,5,6,7]

-- 4. Napisz funkcje
elimMaybe :: c -> (a -> c) -> Maybe a -> c
elimMaybe d m v = case v of { Nothing -> d; Just x -> m x }
fromMaybe :: a -> Maybe a -> a
fromMaybe d v = case v of { Nothing -> d; Just x -> x }
mapMaybe :: (a -> b) -> Maybe a -> Maybe b
mapMaybe m v = case v of { Nothing -> Nothing; Just x -> Just $ m x }
maybeHead :: [a] -> Maybe a
maybeHead t = case t of { [] -> Nothing; (x:_) -> Just x }
elimEither :: (a  -> c) -> (b -> c) -> Either a b -> c
elimEither am bm e = case e of { Left a -> am a; Right b -> bm b}
mapEither :: (a1 -> a2) -> (b1 -> b2) -> Either a1 b1 -> Either a2 b2
mapEither am bm e = case e of { Left a -> Left $ am a; Right b -> Right $ bm b }
mapRight ::  (b1 -> b2) -> Either a b1 -> Either a b2
mapRight m e = case e of { Left x -> Left x; Right x -> Right $ m x }
fromEither :: Either a a -> a
fromEither e = case e of { Left x -> x; Right x -> x }

reverseRight :: Either e [a] -> Either e [a]
reverseRight (Left x) = Left x
reverseRight (Right t) = Right $ reverse t




data Exp
  = EInt Int             -- stała całkowita
  | EAdd Exp Exp         -- e1 + e2
  | ESub Exp Exp         -- e1 - e2
  | EMul Exp Exp         -- e1 * e2
  | EVar String          -- zmienna
  | ELet String Exp Exp  -- let var = e1 in e2

instance Eq Exp where
  (EInt n1) == (EInt n2) = n1 == n2
  (EAdd e1 e2) == (EAdd e1' e2') = e1 == e1' && e2 == e2'
  (ESub e1 e2) == (ESub e1' e2') = e1 == e1' && e2 == e2'
  (EMul e1 e2) == (EMul e1' e2') = e1 == e1' && e2 == e2'
  (EVar v1) == (EVar v2) = v1 == v2
  (ELet v1 e1 e2) == (ELet v2 e1' e2') = v1 == v2 && e1 == e1' && e2 == e2'
  _ == _ = False

instance Show Exp where
  show (EInt n) = show n
  show (EAdd e1 e2) = "(" ++ show e1 ++ " + " ++ show e2 ++ ")"
  show (ESub e1 e2) = "(" ++ show e1 ++ " - " ++ show e2 ++ ")"
  show (EMul e1 e2) = "(" ++ show e1 ++ " * " ++ show e2 ++ ")"
  show (EVar v) = v
  show (ELet v e1 e2) = "let " ++ v ++ " = " ++ show e1 ++ " in " ++ show e2

instance Num Exp where
  (+) = EAdd
  (-) = ESub
  (*) = EMul
  negate e = 0 - e
  abs _ = undefined
  signum _ = undefined
  fromInteger n = EInt (fromInteger n)

testExp2 :: Exp
testExp2 = (2 + 2) * 3

simpl :: Exp -> Exp
simpl (EAdd e1 e2) = let e1' = simpl e1
                         e2' = simpl e2
                     in case (e1', e2') of
                          (EInt 0, _) -> e2'
                          (_, EInt 0) -> e1'
                          _           -> EAdd e1' e2'
simpl (ESub e1 e2) = let e1' = simpl e1
                         e2' = simpl e2
                     in case (e1', e2') of
                          (_, EInt 0) -> e1'
                          _           -> ESub e1' e2'
simpl (EMul e1 e2) = let e1' = simpl e1
                         e2' = simpl e2
                     in case (e1', e2') of
                          (EInt 0, _) -> EInt 0





main :: IO ()
main = do
  putStrLn "Hello, Haskell!"
