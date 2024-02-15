module SparsePoly(fromDP, toDP, qrP) where
import PolyClass
import Representation
import Data.Maybe (listToMaybe)

import Data.List (sortBy)
import Data.Function (on)
import Debug.Trace (trace)

type Poly = [(Float,Int)]

mySort :: Ord a => [(a, b)] -> [(a, b)]
mySort = sortBy (flip compare `on` fst)

polyDegree :: [(Int, a)] -> Int
polyDegree = maybe (-1) fst . listToMaybe

convert :: (Eq a,Num a) => Int -> [(Int, a)] -> [a]
convert (-1) [] = []
convert (-1) _ = error "convert: negative exponent"
convert n [] = 0 : convert (n - 1) []
convert n p@((exp, coeff):rest) = 
    if exp == n then coeff : convert (n - 1) rest
    else 0 : convert (n - 1) p

sparseToDense :: (Eq a, Num a) => SparsePoly a -> DensePoly a
sparseToDense (S t) = 
  if null t then P [] 
  else let max_exp = fst $ head t in
    P $ reverse $ convert max_exp $ normalize t


-- | fromDP example
-- >>> fromDP sampleDP
-- S {unS = [(3,1),(0,-1)]}
fromDP :: (Eq a, Num a) => DensePoly a -> SparsePoly a
toDP :: (Eq a, Num a) => SparsePoly a -> DensePoly a

fromDP (P t) = S $ normalize [(exp, coef) | (exp, coef) <- zip [0..] t, coef /= 0]
toDP = sparseToDense

first :: (a -> a') -> (a, b) -> (a', b)
first m (a, b) = (m a, b)
second :: (b -> b') -> (a, b) -> (a, b')
second m (a, b)= (a, m b)

expFMap :: (Int -> Int) -> SparsePoly a -> SparsePoly a
expFMap m = S . map (first m) . unS

accumulate :: (Eq a, Num a) => [(Int, a)] -> [(Int, a)]
accumulate [] = []
accumulate [x] = [x]
accumulate ((exp1, coef1):(exp2, coef2):r)
    | exp1 == exp2 = accumulate $ (exp1, coef1 + coef2) : r
    | otherwise    = (exp1, coef1) : accumulate ((exp2, coef2):r)

normalize :: (Num a, Eq a) => [(Int, a)] -> [(Int, a)]
normalize = accumulate . filter ((/= 0) . snd) . mySort

instance Functor SparsePoly where
  fmap m = S . map (second m) . unS

instance Polynomial SparsePoly where
  zeroP = S []
  constP w = S [(0, w) | w /= 0]
  varP = S [(1, 1)]
  evalP (S t) v = foldr (\ (exp, coef) acc -> acc + v^exp * coef) 0 t
  shiftP n = expFMap (+n)
  degree = polyDegree . unS

addPoly :: (Eq a, Num a) => [(Int, a)] -> [(Int, a)] -> [(Int, a)]
addPoly [] [] = []
addPoly [] t = t
addPoly t [] = t
addPoly ((exp1, coef1):r1) ((exp2, coef2):r2)
    | exp1 == exp2 = let sum = coef1 + coef2 in
         let rest = addPoly r1 r2 in 
            if sum == 0 then rest 
            else (exp1, sum) : rest
    | exp1 > exp2  = (exp1, coef1) : addPoly r1 ((exp2, coef2):r2)
    | otherwise    = (exp2, coef2) : addPoly ((exp1, coef1):r1) r2

(|>) :: (a -> b) -> (b -> c) -> (a -> c)
(|>) = flip (.)

mulPoly :: (Eq a, Num a) => [(Int, a)] -> [(Int, a)] -> [(Int, a)]
mulPoly [] _ = []
mulPoly _ [] = []
mulPoly ((exp1, coef1):t1) p2 = 
    normalize $ map (\ (exp2, coef2) -> (exp2 + exp1, coef2 * coef1)) p2 ++ mulPoly t1 p2

instance (Eq a, Num a) => Num (SparsePoly a) where
  (+) (S t1) (S t2) = S . normalize $ addPoly t1 t2
  (*) (S t1) (S t2) = S . normalize $ mulPoly t1 t2
  abs = undefined
  signum = undefined
  fromInteger = constP . fromInteger
  negate p = constP (-1) * p

instance (Eq a, Num a) => Eq (SparsePoly a) where
    p == q = nullP(p-q)


-- qrP s t | not(nullP t) = (q, r) iff s == q*t + r && degree r < degree t
-- Polynomial long division
qrP :: (Eq a, Fractional a) => SparsePoly a -> SparsePoly a -> (SparsePoly a, SparsePoly a)
qrP n' d'
  | nullP d' = error "qrP: division by zero"
  | otherwise = let (n, d) = (S $ normalize $ unS n', S $ normalize $ unS d') in qrSHelper n d (0, n) where
    qrSHelper :: (Eq a, Fractional a) => SparsePoly a -> SparsePoly a -> (SparsePoly a, SparsePoly a) -> (SparsePoly a, SparsePoly a)
    qrSHelper n d (q, r)
      | nullP r || degree r < degree d = (q, r)
      | otherwise = qrSHelper n d (q + q', r - r')
      where
        leadCoeff = snd . head . unS
        q' = shiftP (degree r - degree d) (constP (leadCoeff r / leadCoeff d))
        r' = d * q'
          
-- | Division example
-- >>> let x = varP in qrP (x^2 - 1) (x -1) == ((x + 1), 0)
-- True
