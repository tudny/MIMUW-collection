module PolyClass where

class Polynomial p where
  zeroP  :: p a                           -- zero polynomial
  constP :: (Eq a, Num a) => a -> p a     -- constant polynomial
  varP   :: Num a => p a                  -- p(x) = x
  x      :: Num a => p a                  -- p(x) = x
  x       = varP
  evalP  :: Num a => p a -> a -> a        -- value of p(x) at given x
  shiftP :: (Eq a, Num a) => Int -> p a -> p a    -- multiply by x^n
  degree :: (Eq a, Num a) => p a -> Int   -- highest power with nonzero coefficient
  nullP  :: (Eq a, Num a) => p a -> Bool  -- True for zero polynomial
  nullP p = degree p < 0
  