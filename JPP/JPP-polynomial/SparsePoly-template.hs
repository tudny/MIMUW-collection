module SparsePoly(fromDP, toDP, qrP) where
import PolyClass
import Representation

-- | fromDP example
-- >>> fromDP sampleDP
-- S {unS = [(3,1),(0,-1)]}
fromDP :: (Eq a, Num a) => DensePoly a -> SparsePoly a
toDP :: (Eq a, Num a) => SparsePoly a -> DensePoly a

fromDP = undefined
toDP = undefined

first :: (a -> a') -> (a, b) -> (a', b)
first = undefined
second :: (b -> b') -> (a, b) -> (a, b')
second = undefined

instance Functor SparsePoly where

instance Polynomial SparsePoly where

instance (Eq a, Num a) => Num (SparsePoly a) where

instance (Eq a, Num a) => Eq (SparsePoly a) where
    p == q = nullP(p-q)

-- qrP s t | not(nullP t) = (q, r) iff s == q*t + r && degree r < degree t
qrP :: (Eq a, Fractional a) => SparsePoly a -> SparsePoly a -> (SparsePoly a, SparsePoly a)
qrP = undefined

-- | Division example
-- >>> qrP (x^2 - 1) (x -1) == ((x + 1), 0)
-- True