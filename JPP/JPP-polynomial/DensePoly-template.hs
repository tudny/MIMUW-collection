module DensePoly() where
import PolyClass
import Representation

instance Functor DensePoly where

instance Polynomial DensePoly where

instance (Eq a, Num a) => Num (DensePoly a) where

-- |
-- >>> x^3 - 1 :: DensePoly Integer 
-- P {unP = [-1,0,0,1]}

-- | Num operations give canonical results:
-- >>> isCanonicalDP (sampleDP - sampleDP)
-- True
    
instance (Eq a, Num a) => Eq (DensePoly a) where

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