{-# LANGUAGE TemplateHaskell #-}
import Test.QuickCheck
import Data.List(sortOn)
import DensePoly
import SparsePoly
import Representation
import PolyClass
import Debug.Trace (trace)

type DPI = DensePoly Int
type SPI = SparsePoly Int

prop_AddCommDP :: DPI -> DPI -> Property
prop_AddCommDP p q = p + q === q + p

prop_AddZeroRDP :: DensePoly Int -> Property
prop_AddZeroRDP p = p + zeroP === p

prop_MulZeroRDP :: DensePoly Int -> Property
prop_MulZeroRDP p = p * zeroP === zeroP

prop_MulZeroLDP :: DensePoly Int -> Property
prop_MulZeroLDP p = zeroP * p === zeroP

prop_MulCommDP :: DPI -> DPI -> Property
prop_MulCommDP p q = p * q === q * p

prop_OneRDP :: DensePoly Int -> Property
prop_OneRDP p = -- label ("Degree of p is " ++ show(degree p)) $
  p * constP 1 === p

prop_DistLDP :: DPI -> DPI -> DPI -> Property
prop_DistLDP p q r = p*(q+r) === p*q + p*r

prop_MulCanonicalDP :: DPI -> DPI -> Bool
prop_MulCanonicalDP p q = isCanonicalDP (p*q)

prop_SubCanonicalDP :: DPI -> DPI -> Bool
prop_SubCanonicalDP p q = isCanonicalDP (p-q)

prop_ShiftLDP :: NonNegative(Small Int) -> DPI -> DPI -> Property
prop_ShiftLDP (NonNegative (Small n)) p q = shiftP n p * q === shiftP n (p*q)

prop_EvalPlus  :: Int ->  DPI -> DPI -> Property
prop_EvalPlus x p q = evalP(p + q) x === evalP p x + evalP q x

-- SPI

prop_AddCommSP :: SPI -> SPI -> Property
prop_AddCommSP p q = within 100000 $ p + q === q + p

prop_AddZeroRSP :: SPI -> Property
prop_AddZeroRSP p = isCanonicalSP p ==>
                    p + zeroP === p

prop_MulZeroRSP :: SPI -> Property
prop_MulZeroRSP p = isCanonicalSP p ==>
                    p * zeroP === zeroP

prop_MulZeroLSP :: SPI -> Property
prop_MulZeroLSP p = isCanonicalSP p ==>
                    zeroP * p === zeroP

prop_OneRSP :: SPI -> Property
prop_OneRSP p = isCanonicalSP p ==>
                p * constP 1 === p

-- within: prop fails if it does not complete within the given number of microseconds.
prop_MulCommSP :: SPI -> SPI -> Property
prop_MulCommSP p q = within 100000 $ p * q === q * p

prop_DistLSP :: SPI -> SPI -> SPI -> Property
prop_DistLSP p q r = within 100000 $ p*(q+r) === p*q + p*r

prop_ShiftLSP :: NonNegative(Small Int) -> SPI -> SPI -> Property
prop_ShiftLSP (NonNegative (Small n)) p q = shiftP n p * q === shiftP n (p*q)

prop_MulCanonicalSP :: SPI -> SPI -> Bool
prop_MulCanonicalSP p q = isCanonicalSP (p*q)

prop_SubCanonicalSP :: SPI -> SPI -> Bool
prop_SubCanonicalSP p q = isCanonicalSP (p-q)

-- conversions

prop_fromToDP :: SPI -> Property
prop_fromToDP p = isCanonicalSP p ==> fromDP(toDP p) == p

prop_toFromDP :: DPI -> Property
prop_toFromDP p = isCanonicalDP p ==> toDP(fromDP p) == p


type SPR = SparsePoly Rational

prop_qr1 :: SPR -> (NonZero SPR) -> Bool
prop_qr1 p (NonZero s) = p == q*s + r where (q,r) = qrP p s

prop_qr2 :: SPR -> (NonZero SPR) -> Bool
prop_qr2 p (NonZero s) = degree r < degree s where (q,r) = qrP p s

writeln :: String -> IO ()
writeln = putStrLn

-- Hic sunt leones

instance (Num a, Arbitrary a) => Arbitrary (DensePoly a) where
  arbitrary = P <$> arbitrary
  shrink = map P . shrink . unP

log2 :: Int -> Int
log2 0 = 0
log2 n = 1 + log2 (div n 2)

instance (Num a, Eq a, Arbitrary a) => Arbitrary (SparsePoly a) where
  arbitrary = S . norm <$> sized g where
    norm = sortOn (negate . fst)
    g 0 = return []
    g n = do
      let p = log2 n `div` 2
      a <- frequency [(n-p, return 0), (p, arbitrary)]
      r <- g(n-1)
      return $ if a /= 0 then (n,a):r else r
  shrink (S ps) = map S $ s ps where
    s [] = []
    s ((a,n):ps) = ps:[(a,n):ps' | S ps' <- shrink (S ps)]



return []
runTests = $quickCheckAll

main = runTests
