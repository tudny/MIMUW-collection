module Main where
import Data.Char


-- Proste zadania
countdown :: Int -> [Int]
countdown n = [n, n-1 .. 0]


collatz :: Int -> [Int]
collatz n
  | n == 1 = [1]
  | even n = n : collatz (n `div` 2)
  | odd n = n : collatz (3 * n + 1)
  | otherwise = error "Number has to be even or odd"

-- >>> countdown 10
-- [10,9,8,7,6,5,4,3,2,1,0]

-- >>> collatz 10
-- [10,5,16,8,4,2,1]


-- Zadania (raczej po wykładzie)
-- 1. Napisz własne odpowiedniki standardowych funkcji head, tail, ++, take, drop, filter, map, concat

myhead :: [a] -> a
myhead (x:_) = x
myhead [] = error "Empty list"

-- >>> myhead [1,2,3]
-- 1

mytail :: [a] -> [a]
mytail (_:xs) = xs
mytail [] = error "Empty list"

-- >>> mytail [1,2,3]
-- [2,3]

(+++) :: [a] -> [a] -> [a]
(+++) [] t = t
(+++) (h:t) s = h : (t +++ s)

-- >>> [1,2,3] +++ [4,5,6]
-- [1,2,3,4,5,6]

mytake :: Int -> [a] -> [a]
mytake 0 _ = []
mytake n (x:xs) = x : mytake (n - 1) xs
mytake n [] = error "Empty list"

-- >>> mytake 3 [1..10]
-- [1,2,3]

-- >>> mytake 3 [1..2]
-- Empty list

-- >>> mytake 5 [1..]
-- [1,2,3,4,5]

mydrop :: Int -> [a] -> [a]
mydrop 0 t = t
mydrop n (_:xs) = mydrop (n - 1) xs
mydrop n [] = error "Empry list"

-- >>> mydrop 3 [1..10]
-- [4,5,6,7,8,9,10]

myfilter :: (a -> Bool) -> [a] -> [a]
myfilter _ [] = []
myfilter p (x:xs) =
    let rest = myfilter p xs
     in if p x
        then x : rest
        else rest

-- >>> myfilter even [1..10]
-- [2,4,6,8,10]

mymap :: (a -> b) -> [a] -> [b]
mymap _ [] = []
mymap m (x:xs) = m x : mymap m xs

-- >>> mymap succ [1..3]
-- [2,3,4]

myconcat :: [[a]] -> [a]
myconcat [] = []
myconcat (x:xs) = x +++ myconcat xs

-- >>> myconcat ["abc", "def", "ghj"]
-- "abcdefghj"


-- 2. Napisz funkcję inits, ktora dla danej listy da listę wszystkich jej odcinków początkowych, np.
-- inits [1,2] == [[],[1],[1,2]]
inits :: [a] -> [[a]]
inits [] = [[]]
inits (x:xs) = [] : map (x:) (inits xs)

-- >>> inits [1,2,3]
-- [[],[1],[1,2],[1,2,3]]

-- 3. Napisz funkcje partitions, ktora dla danej listy xs da liste wszystkich par (ys,zs) takich, że
-- xs == ys ++ zs
partitions :: [a] -> [([a], [a])]
partitions [] = [([], [])]
partitions (x:xs) = ([], x:xs) : map (\(ys, zs) -> (x:ys, zs)) (partitions xs)

-- >>> partitions [1,2,3]
-- [([],[1,2,3]),([1],[2,3]),([1,2],[3]),([1,2,3],[])]


-- 4. Napisz funkcje permutations, która dla danej listy da listę wszystkich jej permutacji (dla uniknięcia niejasności możemy założyć, ze wszystkie elementy listy wejściowej są różne)

insert :: a -> ([a], [a]) -> [a]
insert x (t1, t2) = t1 ++ [x] ++ t2

permutations :: [a] -> [[a]]
permutations [] = [[]]
permutations (x:xs) =
        concatMap (map (insert x) . partitions) (permutations xs)

-- >>> permutations [1,2,3]
-- [[1,2,3],[2,1,3],[2,3,1],[1,3,2],[3,1,2],[3,2,1]]



-- 5. Napisz funkcje nub, ktora usunie z listy wszystkie duplikaty, np nub [1,2,1,3,1,2,1,4] == [1,2,3,4]
nub :: Eq a => [a] -> [a]
nub [] = []
nub (x:xs) = x : nub (filter (x/=) xs)

-- >>> nub [1,2,1,3,1,2,1,4] == [1,2,3,4]
-- >>> nub [1,2,1,3,1,2,1,4]
-- True
-- [1,2,3,4]



-- 1*. Ćwiczymy funkcje wyzszego rzedu
--       a. Napisz funkcje

-- incAll :: [[Int]] -> [[Int]]
-- która zwiększy o 1 każdy element każdego elementu swojego argumentu, np

-- *Main Data.List> incAll $ inits [1..3] [[],[2],[2,3],[2,3,4]]

incAll :: [[Int]] -> [[Int]]
incAll = map (map (+1))

-- >>> incAll $ inits [1..3]
-- [[],[2],[2,3],[2,3,4]]


--      b. Napisz przy pomocy foldr
--      silnię
--      concat :: [[a]] -> [a]

factorial :: Int -> Int
factorial x = foldr (*) 1 [1..x]

mysecondconcat :: [[a]] -> [a]
mysecondconcat = foldr (++) []

-- >>> factorial 4
-- 24

-- >>> mysecondconcat ["abc", "def", "ghj"]
-- "abcdefghj"


-- 1d. d. Napisz funkcję obliczającą iloczyn skalarny dwóch list liczb; użyj zipWith

iloczynSkalary :: [Int] -> [Int] -> [Int]
iloczynSkalary = zipWith (*)

-- 2. Napisz funkcję triples :: Int -> [(Int,Int,Int)], która dla argumentu n da listę wszystkich trójek liczb laturalnych o elementach z [1..n]

triples :: Int -> [(Int, Int, Int)]
triples n = let seg = [1..n] in [(x, y, z) | x <- seg, y <- seg, z <- seg]

-- >>> triples 2
-- [(1,1,1),(1,1,2),(1,2,1),(1,2,2),(2,1,1),(2,1,2),(2,2,1),(2,2,2)]

-- 3. Napisz funkcję triads :: Int -> [(Int,Int,Int)], ktora da liste trojek pitagorejskich (𝑥2+𝑦2=𝑧2 dla x,y,z <= n)

(|>) :: (a -> b) -> (b -> c) -> a -> c
(|>) = flip (.)

gcdAll :: Integral a => [a] -> a
gcdAll = foldr gcd 0 

triads :: Int -> [(Int, Int, Int)]
triads = filter (\(x, y, z) -> x*x + y*y == z*z) . triples

-- 3* Nietrywialne trójki pitagorejskie
nonTrivialTriads :: Int -> [(Int, Int, Int)]
-- nonTrivialTriads = triads |> filter (\ (x, y, z) -> gcdAll [x, y, z] == 1 && x < y)
nonTrivialTriads n = let sq = [1..n] in
    [(x, y, z) | x <- sq, y <- sq, z <- sq, x*x + y*y == z*z, x < y, gcdAll [x, y, z] == 1]

-- >>> triads 20
-- [(3,4,5),(4,3,5),(5,12,13),(6,8,10),(8,6,10),(8,15,17),(9,12,15),(12,5,13),(12,9,15),(12,16,20),(15,8,17),(16,12,20)]

-- >>> nonTrivialTriads 20
-- [(3,4,5),(5,12,13),(8,15,17)]

-- 4. a. Napisz funkcję incMaybe :: Maybe Int -> Maybe Int

incMaybe :: Maybe Int -> Maybe Int
-- incMaybe = fmap (+1)
incMaybe x = 
    case x of
        Nothing -> Nothing
        Just v -> Just (v + 1)

-- >>> incMaybe (Just 41)
-- Just 42

-- >>> incMaybe Nothing
-- Nothing

-- 4. b. Napisz funkcję addMaybe :: Maybe Int -> Maybe Int -> Maybe Int

opMaybe :: (a -> b -> c) -> (Maybe a -> Maybe b -> Maybe c)
-- opMaybe op x y = fmap op x <*> y
opMaybe op x y = case (x, y) of
    (Just a, Just b) -> Just $ op a b
    (_, _)           -> Nothing

addMaybe :: Maybe Int -> Maybe Int -> Maybe Int
addMaybe = opMaybe (+)

-- >>> addMaybe (Just 4) (Just 6)
-- >>> addMaybe Nothing (Just 6)
-- >>> addMaybe (Just 6) Nothing
-- >>> addMaybe Nothing Nothing
-- Just 10
-- Nothing
-- Nothing
-- Nothing

-- 5. E[| skip |]
-- 5. Można zrobic

-- Liczby Fibonacciego
-- lista wszystkich liczb Fibonacciego
fibs = 1 : 1 : zipWith (+) fibs (tail fibs)

-- liczenie liczb pierwszych na rozne sposoby
primes = sive [2..]
    where 
        sive [] = []
        sive (x:xs) = x : [y | y <- xs, y `mod` x /= 0]

-- >>> take 10 primes
-- [2,3,5,7,9,11,13,15,17,19]

-- Rekursja ogonowa -- silnia z akumulatorem, liczby Fibbonacciego, reverse
facAcc n = facAccHelper n 1
    where
        facAccHelper 0 acc = acc
        facAccHelper n acc = facAccHelper (n - 1) (acc * n)

-- >>> facAcc 4
-- 24



-- >>> take 10 fibs
-- [1,1,2,3,5,8,13,21,34,55]

-- 6. Napisz funkcję

-- indexOf :: Char -> String -> Maybe Int
-- taką, że jeśli c wystepuje w s na pozycji n, to (indexOf c s) == Just n, wpp Nothing, np

-- *Main> indexOf 'a' "Ala"
-- Just 2
-- *Main> indexOf 'b' "Ala"
-- Nothing

indexOf :: Char -> String -> Maybe Int
indexOf c s = indexOfHelper c s 0
    where
        indexOfHelper :: Char -> String -> Int -> Maybe Int
        indexOfHelper c [] _ = Nothing
        indexOfHelper c (s:ss) id = if c == s then Just id else indexOfHelper c ss (id + 1)

-- >>> indexOf 'a' "Ala"
-- >>> indexOf 'b' "Ala"
-- Just 2
-- Nothing

-- Napisz funkcje

-- positions :: Char -> String -> [Int]
-- taką, że (positions c s) daje listę wszystkich pozycji, na których c wystepuje w s, np.

-- *Main> positions 'a' "Ala ma kota"
-- [2,5,10]
-- *Main> positions 'b' "Ala ma kota"
-- []

positions :: Char -> String -> [Int]
positions c s = fst $ foldr (\cc (r, id) -> if c == cc then (id:r, id - 1) else (r, id - 1)) ([], length s - 1) s

positions2 :: Char -> String -> [Int]
positions2 c s = positions2Helper c s 0
    where
        positions2Helper :: Char -> String -> Int -> [Int]
        positions2Helper c [] _ = []
        positions2Helper c (s:ss) id = let rest = positions2Helper c ss (id + 1) in
            if s == c then id:rest else rest

-- >>> positions 'a' "Ala ma kota"
-- >>> positions 'b' "Ala ma kota"
-- [2,5,10]
-- []

-- >>> positions2 'a' "Ala ma kota"
-- >>> positions2 'b' "Ala ma kota"
-- [2,5,10]
-- []

showInt :: Int -> String
showInt 0 = "0"
showInt n =
    let dgToChr d = chr (d + ord '0') in
    let digit = [dgToChr (n `mod` 10)] in
    let rest = n `div` 10 in
        if rest /= 0 then
            showInt rest ++ digit
        else
            digit

-- >>> showInt 123
-- "123"

-- >>> showInt 0
-- "0"

showLst :: (a -> String) -> [a] -> String
showLst m l = "[" ++ drop 2 (concatMap (\x -> ", " ++ m x) l) ++ "]"

showIntLst :: [Int] -> String
showIntLst = showLst showInt

-- >>> showIntLst [123, 54, 1, 0, 99, 3]
-- "[123, 54, 1, 0, 99, 3]"

-- 8. Napisz program, który podzieli tekst otrzymany na wejściu na linie długości nie większej niż zadana stała (np. 30)

splitOf :: Int -> String -> [String]
splitOf c [] = []
splitOf c s = let (l, r) = splitAt c s in l : splitOf c r

-- >>> splitOf 3 "Ala ma kota, kot ma Ale"
-- ["Ala"," ma"," ko","ta,"," ko","t m","a A","le"]

-- >>> splitOf 4 "a"
-- ["a"]

-- Main

main :: IO ()
main = do
    print (countdown 10)
    print (collatz 10)


myflip :: (a -> b -> c) -> (b -> a -> c)
myflip f x y = f y x

colon = myflip (.)

f = (+2)
g = (*3)

-- >>> f . g $ 5
-- >>> f `colon` g $ 5
-- 17
-- 21

