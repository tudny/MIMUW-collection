module ReadInts where
import Data.Char (isDigit)
-- 6. Stwórz moduł ReadInts
-- a. Napisz funkcję

-- readInts :: String -> [Int]
-- która odczyta z napisu występujące w nim liczby naturalne, np

-- *Main> readInts "1 23 456 7.8 abc 9"
-- [1,23,456,9]
-- *Main> readInts "foo"
-- []
-- użyj funkcji isDigit z modulu Data.Char oraz funkcji map, filter, all z Prelude

-- b. Napisz podobną funkcję readInts2 :: String -> Either String [Int] która da listę liczb, jeśli wszystkie słowa jej argumentu są liczbami a komunikat o błędzie w przeciwnym przypadku

-- Może się przydać funkcja reverseRight (albo mapRight)

--     *Main> readInts2  "1 23 456 foo 9"
--     Left "Not a number: foo"
--     *Main> readInts2  "1 23 456"     
--     Right [1,23,456]
-- c. Napisz funkcję

-- sumInts :: String -> String
-- jesli wszystkie slowa jej argumentu są liczbami da reprezentacje ich sumy
-- wpp komunikat o bledzie
-- stwórz program importujący stworzony moduł i uruchamiający funkcję sumInts przy pomocy interact.


-- Global helpers

(|>) :: (a -> b) -> (b -> c) -> a -> c
(|>) = flip (.)

mapRight ::  (b1 -> b2) -> Either a b1 -> Either a b2
mapRight m e = case e of { Left x -> Left x; Right x -> Right $ m x }

readEither :: String -> Either String Int
readEither s = if all isDigit s then Right (read s) else Left ("Not a number: " ++ s)


-- a)
readInts :: String -> [Int]
readInts = words |> filter (all isDigit) |> map read

-- >>> readInts "1 23 456 7.8 abc 9"
-- [1,23,456,9]

-- b)
readInts2 :: String -> Either String [Int]
readInts2 = words 
    |> map readEither 
    |> foldr (\x acc -> case (x, acc) of 
        (Right x, Right acc) -> Right (x:acc)
        (Left x, _) -> Left x
        (_, Left acc) -> Left acc) (Right [])

-- >>> readInts2  "1 23 456 foo 9"
-- Left "Not a number: foo"

-- >>> readInts2  "1 23 456"
-- Right [1,23,456]

