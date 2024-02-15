module Representation where

-- | Polynomial as a list of coefficients, lowest power first
-- e.g. x^3-1 is represented as P [-1,0,0,1]
-- canonical: no trailing zeros
newtype DensePoly a = P { unP :: [a] } deriving Show

isCanonicalDP :: (Eq a, Num a) => DensePoly a -> Bool
isCanonicalDP = go . unP where
    go [] = True
    go xs = last xs /= 0

sampleDP = P [-1,0,0,1]

-- | Polynomial as a list of coefficients, highest power first
-- e.g. x^3-1 is represented as R [1,0,0,-1]
-- canonical: no leading zeros
newtype ReversePoly a = R { unR :: [a] } deriving Show

isCanonicalRP :: (Eq a, Num a) => ReversePoly a -> Bool
isCanonicalRP = go . unR where
    go [] = True
    go (x:_) = x /= 0
    
sampleRP = R [1,0,0,-1]

-- | Polynomial as a list of pairs (power, coefficient)
-- e.g. x^3-1 is represented as S [(3,1),(0,-1)]
-- canonical: in descending order of powers; no zero coefficients
newtype SparsePoly a = S { unS :: [(Int, a)] } deriving Show


decreasing :: Ord a => [a] -> Bool
decreasing [] = True
decreasing [_] = True
decreasing (x:xs@(y:ys)) = x > y && decreasing xs
                           
isCanonicalSP (S ps) = decreasing (map fst ps) && all ((/=0).snd) ps

sampleSP = S [(3,1),(0,-1)]