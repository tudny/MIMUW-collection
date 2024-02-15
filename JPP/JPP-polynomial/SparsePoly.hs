module SparsePoly(fromDP, toDP, qrP) where
import PolyClass
import Representation
import Data.List (dropWhileEnd)
import Data.Function (on)
import Data.Maybe (listToMaybe, fromMaybe)

-- | fromDP example
-- >>> fromDP sampleDP
-- S {unS = [(3,1),(0,-1)]}
fromDP :: (Eq a, Num a) => DensePoly a -> SparsePoly a
toDP :: (Eq a, Num a) => SparsePoly a -> DensePoly a

fromDP = S . normalizeInternal . denseToSparseInternal . normalizeInternalDense . unP
toDP = P . normalizeInternalDense . sparseToDenseInternal . normalizeInternal . unS

first :: (a -> a') -> (a, b) -> (a', b)
first f (a, b) = (f a, b)
second :: (b -> b') -> (a, b) -> (a, b')
second f (a, b) = (a, f b)

instance Functor SparsePoly where
  fmap = wrapped . map . second

instance Polynomial SparsePoly where
  zeroP = S polyZero
  constP = S . polyConst
  varP = S polyLinear
  evalP = polyEval . unS
  shiftP = wrapped . polyShift
  degree = polyDegree . unS

instance (Eq a, Num a) => Num (SparsePoly a) where
  (+) = normalizeBiOp polyAdd
  (*) = normalizeBiOp polyMul
  abs = undefined
  signum = undefined
  fromInteger = constP . fromInteger
  negate = normalizeOp polyNeg


instance (Eq a, Num a) => Eq (SparsePoly a) where
    p == q = nullP (p - q)

-- qrP s t | not(nullP t) = (q, r) iff s == q*t + r && degree r < degree t
qrP :: (Eq a, Fractional a) => SparsePoly a -> SparsePoly a -> (SparsePoly a, SparsePoly a)
qrP n' d'
  | nullP d = error "qrP: division by zero"
  | otherwise = qrSHelper (0, n)
  where
    n = normalize n'
    d = normalize d'
    leadCoeff = snd . head . unS
    qrSHelper (q, r)
      | nullP r || degree r < degree d = (q, r)
      | otherwise = qrSHelper (q + q', r - r')
      where
        q' = shiftP (degree r - degree d) (constP (leadCoeff r / leadCoeff d))
        r' = d * q'

-- | Division example
-- >>> qrP (x^2 - 1) (x -1) == ((x + 1), 0)
-- True

-- ============================================================================
-- | Sort

sortByInternal :: (a -> a -> Ordering) -> [a] -> [a]
sortByInternal _ [] = []
sortByInternal cmp (x:xs) = left ++ midd ++ right
  where
    left = sortByInternal cmp [y | y <- xs, cmp y x == LT]
    midd = x : [y | y <- xs, cmp y x == EQ]
    right = sortByInternal cmp [y | y <- xs, cmp y x == GT]

-- ============================================================================
-- Dense normalization (borrowed from DensePoly.hs)

-- Remove trailing zeros from a list of coefficients
-- foldr is used to avoid traversing the list twice
removeTrailingZerosDense :: (Eq a, Num a) => [a] -> [a]
removeTrailingZerosDense = reverse . dropWhile (== 0) . reverse

-- Normalize a list of coefficients by removing trailing zeros
normalizeInternalDense :: (Eq a, Num a) => [a] -> [a]
normalizeInternalDense = removeTrailingZerosDense

-- ============================================================================
-- Utility functions

-- Sort a polynomial by exponent
sortByExpInternal :: [(Int, a)] -> [(Int, a)]
sortByExpInternal = sortByInternal (flip compare `on` fst)

-- Accumulate similar terms
accumulateInternal :: (Eq a, Num a) => [(Int, a)] -> [(Int, a)]
accumulateInternal [] = []
accumulateInternal [x] = [x]
accumulateInternal ((exp1, coef1):(exp2, coef2):r)
    | exp1 == exp2 = accumulateInternal $ (exp1, coef1 + coef2) : r
    | otherwise    = (exp1, coef1) : accumulateInternal ((exp2, coef2):r)

-- Remove zero terms
removeZeroInternal :: (Eq a, Num a) => [(Int, a)] -> [(Int, a)]
removeZeroInternal = filter ((/= 0) . snd)

-- Normalize a polynomial
normalize :: (Eq a, Num a) => SparsePoly a -> SparsePoly a
normalize = wrapped normalizeInternal

-- We remove zeros twice because first run may remove some terms and we need to remove zeros after accumulation
normalizeInternal :: (Eq a, Num a) => [(Int, a)] -> [(Int, a)]
normalizeInternal = removeZeroInternal . accumulateInternal . sortByExpInternal . removeZeroInternal

-- Run with normalization
normalizeBiOp :: (Eq a, Num a) => ([(Int, a)] -> [(Int, a)] -> [(Int, a)]) -> SparsePoly a -> SparsePoly a -> SparsePoly a
normalizeBiOp op (S p) (S q) = S $ normalizeInternal $ op (normalizeInternal p) (normalizeInternal q)

normalizeOp :: (Eq a, Num a) => ([(Int, a)] -> [(Int, a)]) -> SparsePoly a -> SparsePoly a
normalizeOp op = wrapped $ normalizeInternal . op . normalizeInternal

-- Wrapped operation
wrapped :: ([(Int, a)] -> [(Int, b)]) -> SparsePoly a -> SparsePoly b
wrapped f = S . f . unS

-- ============================================================================
-- Conversion functions

-- Convert a dense polynomial to a sparse polynomial
denseToSparseInternal :: (Eq a, Num a) => [a] -> [(Int, a)]
denseToSparseInternal t = reverse [term | term@(_, coef) <- zip [0..] t, coef /= 0]

-- Convert a sparse polynomial to a dense polynomial
sparseToDenseInternal :: (Eq a, Num a) => [(Int, a)] -> [a]
sparseToDenseInternal = reverse . fst . foldr (\ (exp, coef) (res, len) -> let diff = exp - len in (coef : replicate diff 0 ++ res, len + diff + 1)) ([], 0)

-- ============================================================================
-- Polynomial operations

polyAdd :: (Eq a, Num a) => [(Int, a)] -> [(Int, a)] -> [(Int, a)]
polyAdd [] t = t
polyAdd t [] = t
polyAdd all1@(term1@(exp1, coef1):xs1) all2@(term2@(exp2, coef2):xs2)
  | exp1 == exp2 = append (exp1, sum) $ polyAdd xs1 xs2
  | exp1 > exp2  = term1 : polyAdd xs1 all2
  | otherwise    = term2 : polyAdd all1 xs2
  where
    sum = coef1 + coef2
    append term@(exp, coef) t = [term | coef /= 0] ++ t

polyMul :: (Eq a, Num a) => [(Int, a)] -> [(Int, a)] -> [(Int, a)]
polyMul [] _ = []
polyMul ((exp, coef):rest) q = polyAdd (polyShift exp $ map (second (*coef)) q) (polyMul rest q)

polyNeg :: (Eq a, Num a) => [(Int, a)] -> [(Int, a)]
polyNeg = map (second negate)

polyZero :: [(Int, a)]
polyZero = []

polyConst :: (Eq a, Num a) => a -> [(Int, a)]
polyConst c = [(0, c) | c /= 0]

polyLinear :: Num a => [(Int, a)]
polyLinear = [(1, 1)]

polyEval :: Num a => [(Int, a)] -> a -> a
polyEval p x = sum [coef * x ^ exp | (exp, coef) <- p]

polyShift :: (Eq a, Num a) => Int -> [(Int, a)] -> [(Int, a)]
polyShift = map . first . (+)

polyDegree :: (Eq a, Num a) => [(Int, a)] -> Int
polyDegree = maybe (-1) fst . listToMaybe . normalizeInternal
