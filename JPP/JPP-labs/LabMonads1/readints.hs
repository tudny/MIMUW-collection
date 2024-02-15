import Data.Char (isDigit)
import Text.Read (readEither)
import Control.Monad (foldM, join)
import Control.Applicative (Applicative(liftA2))


fromEither :: Either a a -> a
fromEither = either id id


readInts :: String -> [Int]
readInts = map read . filter (all isDigit) . words


readInts2 :: String -> Either String [Int]
readInts2 = traverse readEither . words


sumInts :: String -> String
sumInts = fromEither . fmap (show . sum) . readInts2


-- >>> readInts2 "1 23 456 abc 9"
-- Left "Prelude.read: no parse"


-- >>> sumInts "1 34 75"
-- >>> sumInts "1 34 abc 75"
-- "110"
-- "Prelude.read: no parse"


-- Zadanie 2

data Exp = Val Int | Div Exp Exp | Add Exp Exp | Sub Exp Exp | Mul Exp Exp
safediv :: Int -> Int -> Maybe Int
safediv _ 0 = Nothing
safediv x y = Just (div x y)

eval :: Exp -> Maybe Int
eval (Val n) = pure n
eval (Div x y) = join $ safediv <$> eval x <*> eval y
eval (Add x y) = (+) <$> eval x <*> eval y
eval (Sub x y) = (-) <$> eval x <*> eval y
eval (Mul x y) = (*) <$> eval x <*> eval y


sampleExp = Div (Div (Val 1972) (Val 2)) (Val 23)
sampleZeroExp = Div (Div (Val 1) (Val 0)) (Val 7)

-- >>> eval sampleExp
-- >>> eval sampleZeroExp
-- Just 42
-- Nothing

evalList' :: [Exp] -> [Maybe Int]
evalList' = map eval

evalList :: [Exp] -> Maybe [Int]
evalList = traverse eval


safediv2 :: Int -> Int -> Either String Int
safediv2 _ 0 = Left "Division by zero"
safediv2 x y = Right (div x y)

eval2 :: Exp -> Either String Int
eval2 (Val n) = pure n
eval2 (Div x y) = join $ safediv2 <$> eval2 x <*> eval2 y
eval2 (Add x y) = (+) <$> eval2 x <*> eval2 y
eval2 (Sub x y) = (-) <$> eval2 x <*> eval2 y
eval2 (Mul x y) = (*) <$> eval2 x <*> eval2 y

