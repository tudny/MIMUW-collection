module DensePoly() where
import PolyClass
import Representation

instance Functor DensePoly where
    fmap m P { unP = x } = P { unP = map m x }

-- remove trailing zeros
removeTrailingZeros :: (Eq a, Num a) => [a] -> [a]
removeTrailingZeros = reverse . dropWhile (== 0) . reverse

normalize :: (Eq a, Num a) => DensePoly a -> DensePoly a
normalize = P . removeTrailingZeros . unP

(|>) = flip (.)

instance Polynomial DensePoly where
  zeroP = P { unP = [] }
  constP x = P { unP = [x | x /= 0] }
  varP = P { unP = [0, 1] }
  evalP P { unP = t } x = case t of
    [] -> 0
    (w:ws) -> w + x * evalP P { unP = ws } x
  shiftP n P { unP = t } = P { unP = replicate n 0 ++ t }
  degree P { unP = t } = length t - 1

zipWithDefault :: (t1 -> t2 -> a) -> t1 -> t2 -> [t1] -> [t2] -> [a]
zipWithDefault op aDef bDef aCol bCol = case (aCol, bCol) of
        (x:xs, y:ys) -> op x y : zipWithDefault op aDef bDef xs ys
        ([], y:ys)   -> op aDef y : zipWithDefault op aDef bDef [] ys
        (x:xs, [])   -> op x bDef : zipWithDefault op aDef bDef xs []
        ([], [])     -> []

polyMul :: (Num a, Eq a) => [a] -> [a] -> [a]
polyMul [] _ = []
polyMul (x:xs) t = zipWithDefault (+) 0 0 (map (*x) t) (unP $ shiftP 1 $ P $ polyMul xs t)

instance (Eq a, Num a) => Num (DensePoly a) where
  (+) P { unP = p } P { unP = q } = normalize $ P { unP = zipWithDefault (+) 0 0 p q }
  (*) (P p) (P q) = normalize $ P $ polyMul p q
  abs = undefined
  signum = undefined
  fromInteger = constP . fromInteger
  negate p = constP (-1) * p

-- |
-- >>> let x = varP :: DensePoly Integer in x^3 - 1
-- P {unP = [-1,0,0,1]}
instance (Eq a, Num a) => Eq (DensePoly a) where
  p == q = unP (normalize p) == unP (normalize q)
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

-- |
-- >>> x^3 - 1 :: DensePoly Integer 
-- P {unP = [-1,0,0,1]}

-- | Num operations give canonical results:
-- >>> isCanonicalDP (sampleDP - sampleDP)
-- True


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
