module DensePoly() where
import PolyClass
import Representation
import Data.List (dropWhileEnd)

instance Functor DensePoly where
  fmap = wrapped . map

instance Polynomial DensePoly where
  zeroP = P polyZero
  constP = P . polyConst
  varP = P polyLinear
  evalP = polyEval . unP
  shiftP n =  normalizeOp $ polyShift n
  degree = polyDegree . normalizeInternal . unP

instance (Eq a, Num a) => Num (DensePoly a) where
  (+) = normalizeBiOp polyAdd
  (*) = normalizeBiOp polyMul
  abs = undefined
  signum = undefined
  fromInteger = constP . fromInteger
  negate = normalizeOp polyNeg

-- |
-- >>> x^3 - 1 :: DensePoly Integer 
-- P {unP = [-1,0,0,1]}

-- | Num operations give canonical results:
-- >>> isCanonicalDP (sampleDP - sampleDP)
-- True

instance (Eq a, Num a) => Eq (DensePoly a) where
    p == q = nullP (p - q)
-- |
-- >>>  P [1,2] == P [1,2]
-- True

-- |
-- >>> fromInteger 0 == (zeroP :: DensePoly Int)
-- True

-- |
-- >>>  P [0,1] == P [1,0]
-- False

-- | Degree examples
-- >>> degree (zeroP :: DensePoly Int)
-- -1
-- >>> degree (constP 1 :: DensePoly Int)
-- 0



-- ============================================================================
-- Utility functions

-- I really liked this operator in OCaml, so I'm bringing it here
infixl 0 |>
(|>) :: a -> (a -> b) -> b
(|>) = flip ($)

-- Zips two lists with a default value for the shorter list
zipWithDefault :: (t1 -> t2 -> a) -> t1 -> t2 -> [t1] -> [t2] -> [a]
zipWithDefault op aDef bDef aCol bCol = case (aCol, bCol) of
        (x:xs, y:ys) -> op x y : zipWithDefault op aDef bDef xs ys
        ([], y:ys)   -> op aDef y : zipWithDefault op aDef bDef [] ys
        (x:xs, [])   -> op x bDef : zipWithDefault op aDef bDef xs []
        ([], [])     -> []

-- ============================================================================
-- Polynomial utility functions

-- Remove trailing zeros from a list of coefficients
-- foldr is used to avoid traversing the list twice
removeTrailingZeros :: (Eq a, Num a) => [a] -> [a]
removeTrailingZeros = dropWhileEnd (== 0)

-- Normalize a list of coefficients by removing trailing zeros
normalizeInternal :: (Eq a, Num a) => [a] -> [a]
normalizeInternal = removeTrailingZeros

-- Normalize a polynomial by removing trailing zeros
normalize :: (Eq a, Num a) => DensePoly a -> DensePoly a
normalize = normalizeInternal |> wrapped

-- Run with normalization
normalizeBiOp :: (Eq a, Num a) => ([a] -> [a] -> [a]) -> DensePoly a -> DensePoly a -> DensePoly a
normalizeBiOp op (P p) (P q) = op (normalizeInternal p) (normalizeInternal q) |> normalizeInternal |> P

normalizeOp :: (Eq a, Num a) => ([a] -> [a]) -> DensePoly a -> DensePoly a
normalizeOp op (P p) = op (normalizeInternal p) |> normalizeInternal |> P

fma :: (Num a) => a -> a -> a -> a
fma x w a = w + a * x

wrapped :: ([a] -> [b]) -> DensePoly a -> DensePoly b
wrapped f = P . f . unP

-- ============================================================================
-- Polynomial operations - all operations operate on the list of coefficients

polyAdd :: (Eq a, Num a) => [a] -> [a] -> [a]
polyAdd = zipWithDefault (+) 0 0

polyMul :: (Num a, Eq a) => [a] -> [a] -> [a]
polyMul [] _ = []
polyMul (x:xs) t = polyAdd (map (*x) t) (polyShift 1 $ polyMul xs t)

polyNeg :: (Num a, Eq a) => [a] -> [a]
polyNeg = map negate

polyZero :: [a]
polyZero = []

polyConst :: (Num a, Eq a) => a -> [a]
polyConst n = [n | n /= 0]

polyLinear :: Num a => [a]
polyLinear = [0, 1]

polyEval :: Num a => [a] -> a -> a
polyEval t x = foldr (fma x) 0 t

polyShift :: (Num a, Eq a) => Int -> [a] -> [a]
polyShift n = (++) (replicate n 0)

polyDegree :: (Num a, Eq a) => [a] -> Int
polyDegree = (+) (-1) . length
