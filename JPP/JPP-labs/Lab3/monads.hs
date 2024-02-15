import Data.Function (on)
import Control.Monad.Reader
import Control.Applicative (Applicative(liftA2))
import Data.Map (Map)
import qualified Data.Map as Map
import Control.Monad.State

-- Zadanie 1

-- Funkcja
allPairs :: [a] -> [a] -> [[a]]
allPairs xs ys = [[x,y] | x <- xs, y <- ys]
-- daje listę wszystkich list dwuelementowych, gdzie pierwszy element należy do pierwszego argumentu, drugi do drugiego, np.

-- >>> allPairs [1,2,3] [4,5]
-- [[1,4],[1,5],[2,4],[2,5],[3,4],[3,5]]


-- przepisz na notację monadyczną (do)
allPairsMonadic :: [a] -> [a] -> [[a]]
allPairsMonadic xs ys = do
    x <- xs
    y <- ys
    [[x, y]]

-- >>> allPairsMonadic [1,2,3] [4,5]
-- [[1,4],[1,5],[2,4],[2,5],[3,4],[3,5]]

-- uogólnij te funkcję do allCombinations :: [[a]] -> [[a]], 
-- która dla n-elementowej listy list da listę wszystkich n-elementowych list takich, 
-- że i-ty element należy do i-tego elementu argumentu, np

-- Main> allCombinations [[1,2], [4,5], [6], [7]]  
-- [[1,4,6,7],[1,5,6,7],[2,4,6,7],[2,5,6,7]]

allCombinations :: [[a]] -> [[a]]
allCombinations = foldr (liftA2 (:)) [[]]

-- >>> allCombinations [[1,2], [4,5], [6], [7]]
-- [[1,4,6,7],[1,5,6,7],[2,4,6,7],[2,5,6,7]]



-- Zadanie 2

-- a.

data Tree a = Empty | Node a (Tree a) (Tree a) deriving (Eq, Ord, Show)
-- Napisz funkcję

renumber :: Tree a -> Tree Int
-- która dla danego drzewa da drzewo, gdzie w każdym węźle przechowywana będzie głębokość tego węzła (odległość od korzenia).
-- Porównaj rozwiązania z użyciem monady Reader i bez.

renumber = _renumber 0
    where
        _renumber :: Int -> Tree a -> Tree Int
        _renumber _ Empty = Empty
        _renumber d (Node _ l r) = (Node d `on` _renumber (d + 1)) l r


-- Using Reader monad
renumberMonadic :: Tree a -> Tree Int
renumberMonadic t = runReader (go t) 0
    where
        go :: Tree a -> Reader Int (Tree Int)
        go Empty = do
            return Empty
        go (Node x l r) = do
            d <- ask
            l' <- local (+1) $ go l
            r' <- local (+1) $ go r
            return $ Node d l' r'

sampleTree = Node 1 (Node 2 Empty Empty) (Node 3 (Node 7 Empty Empty) Empty)
-- >>> renumberMonadic sampleTree
-- Node 0 (Node 1 Empty Empty) (Node 1 (Node 2 Empty Empty) Empty)



-- b. Dane typy dla wyrażeń arytmetycznych

type Var = String
data Exp = EInt Int
     | EOp  Op Exp Exp
     | EVar Var
     | ELet Var Exp Exp  -- let var = e1 in e2

data Op = OpAdd | OpMul | OpSub

type Env = Map Var Int

evalExp :: Exp -> Int
evalExp e = runReader (go e) Map.empty
    where
        run :: Op -> Int -> Int -> Int
        run op = case op of
            OpAdd -> (+)
            OpMul -> (*)
            OpSub -> (-)
        go :: Exp -> Reader Env Int
        go (EInt n) = do
            pure n
        go (EOp op e1 e2) = do
            n1 <- go e1
            n2 <- go e2
            pure $ run op n1 n2
        go (EVar v) = do
            asks (Map.findWithDefault 0 v)
        go (ELet v e1 e2) = do
            n <- go e1
            local (Map.insert v n) $ go e2


sampleExp = ELet "x" (EOp OpAdd (EInt 1) (EInt 2)) (EOp OpMul (EVar "x") (EInt 3))

-- >>> evalExp sampleExp
-- 9




-- Zadanie 3. Monada State

-- a. Dany typ drzew

-- j.w.
-- data Tree a = Empty | Node a (Tree a) (Tree a) deriving (Eq, Ord, Show)

-- Napisz funkcję

renumberTree :: Tree a -> Tree Int
-- która ponumeruje wezly drzewa tak, ze kazdy z nich bedzie mial inny numer. Porownaj rozwiazania z uzyciem monady State i bez.
-- możliwe dodatkowe wymaganie: ponumeruj wezly drzewa w kolejnosci infiksowej.
-- (toList $ renumber $ fromList "Learn Haskell") == [0..12]

-- runState -> (val, env)
-- evalState -> val
renumberTree t = evalState (go t) 0
    where
        go :: Tree a -> State Int (Tree Int)
        go Empty = return Empty
        go (Node x l r) = do
            l' <- go l
            d <- get
            put (d + 1)
            r' <- go r
            return $ Node d l' r'


sampleTree2 = Node 1 (Node 2 Empty Empty) (Node 3 (Node 7 Empty Empty) Empty)

-- >>> renumberTree sampleTree2
-- Node 1 (Node 0 Empty Empty) (Node 3 (Node 2 Empty Empty) Empty)


--                  1
--                 / \
--                0   3
--                   /
--                  2



-- możliwe dodatkowe wymaganie: ponumeruj wezly drzewa w kolejnosci infiksowej.

renumberTreeInfix :: Tree a -> Tree Int
renumberTreeInfix t = evalState (go t) 0
    where
        go :: Tree a -> State Int (Tree Int)
        go Empty = return Empty
        go (Node x l r) = do
            d <- get
            put (d + 1)
            l' <- go l
            r' <- go r
            return $ Node d l' r'

-- >>> renumberTreeInfix sampleTree2
-- Node 0 (Node 1 Empty Empty) (Node 2 (Node 3 Empty Empty) Empty)

--                 0
--                / \
--               1   2
--                  /
--                 3


-- | finite
mergeSort :: (Ord a) => [a] -> [a]
mergeSort t | length t < 2 = t
mergeSort t = merge l r
    where
        mid = length t `div` 2
        enumerated = zip t [0..]
        l = mergeSort [e | (e, n) <- enumerated, n < mid]
        r = mergeSort [e | (e, n) <- enumerated, n >= mid]
        merge :: (Ord a) => [a] -> [a] -> [a]
        merge [] t = t
        merge t [] = t
        merge xx@(x:xs) yy@(y:ys)
            | x < y     = x : merge xs yy
            | otherwise = y : merge xx ys

-- >>> mergeSort [6,56,45,4,568,79,69,675,565,68,6,8,69,9]
-- [4,6,6,8,9,45,56,68,69,69,79,565,568,675]


-- b. Rozszerzmy język z poprzedniego zadania o instrukcje języka Tiny (patrz przedmiot Semantyka i Weryfikacja Programów)

-- Stmt:   S ::= skip | x := e | S1;S2
--         | if b then S1 else S2 | while b do S
-- korzystając z wcześniej napisanej funkcji evalExp, napisz funkcję


data Stmt = SSkip | SAss Var Exp | SCol Stmt Stmt
    | SIf Bool Stmt Stmt | SWhile Bool Stmt

execStmt :: Stmt -> IO ()
-- która wykona podaną instrukcję (program) i wypisze stan końcowy (w tym wypadku wartości zmiennych)


execStmt s = let (_, state) = runState (go s) Map.empty in print state
    where
        go :: Stmt -> State Env ()
        go SSkip = pure ()
        go (SAss v e) = do
            let n = evalExp e
            env <- get
            put (Map.insert v n env)
            pure ()
        go (SCol s1 s2) = do 
            s1' <- go s1
            s2' <- go s2
            pure ()
        go (SIf b s1 s2) 
            | b = do
                s <- go s1
                pure ()
            | otherwise = do
                s <- go s2
                pure ()
        go w@(SWhile b s) 
            | b = do
                s' <- go s
                s' <- go w
                pure ()
            | otherwise = do
                pure ()


sampleStatement = SCol (SAss "x" (EInt 1)) (SAss "y" (EInt 2))


main :: IO ()
main = execStmt sampleStatement

