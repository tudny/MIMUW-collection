import Control.Applicative (Applicative (liftA2))
import Control.Monad (join)
import Data.Function (on)
import qualified Data.Map

data Exp = Val Int | Div Exp Exp | Add Exp Exp | Sub Exp Exp | Mul Exp Exp | Var String | Let String Exp Exp

type Env = Data.Map.Map String Int

safediv :: Int -> Int -> Either String Int
safediv _ 0 = Left "Division by zero"
safediv x y = Right (div x y)

biOperation :: (Int -> Int -> b) -> Env -> Exp -> Exp -> Either String b
-- biOperation op rho x y = op <$> eval rho x <*> eval rho y
-- biOperation op rho = liftA2 op `on` eval rho
biOperation = (. eval) . on . liftA2

maybeToEither :: e -> Maybe a -> Either e a
maybeToEither _ (Just v) = Right v
maybeToEither e Nothing = Left e

eval :: Env -> Exp -> Either String Int
eval rho (Val n) = pure n
eval rho (Div x y) = join $ biOperation safediv rho x y
eval rho (Add x y) = biOperation (+) rho x y
eval rho (Sub x y) = biOperation (-) rho x y
eval rho (Mul x y) = biOperation (*) rho x y
eval rho (Var s) = maybeToEither ("Variable " ++ s ++ " is missing.") (Data.Map.lookup s rho)
eval rho (Let s s_exp e) = do
  s_val <- eval rho s_exp
  eval (Data.Map.insert s s_val rho) e

sampleExprWithAddSubMulDivAndVars =
  Div (Add (Val 1) (Var "x")) (Mul (Sub (Val 10) (Var "y")) (Val 2))

-- >>> eval (Data.Map.fromList [("x", 5), ("y", 3)]) sampleExprWithAddSubMulDivAndVars
-- Right 0

evenBiggerSampleExprWithAddSubMulDivAndVars =
  Div (Add (Val 1) (Var "x")) (Mul (Sub (Val 10) (Var "y")) (Var "z"))

-- >>> eval (Data.Map.fromList [("x", 5), ("y", 3)]) evenBiggerSampleExprWithAddSubMulDivAndVars
-- Left "Variable z is missing."

evalList' :: Env -> [Exp] -> [Either String Int]
evalList' = map . eval

evalList :: Env -> [Exp] -> Either String [Int]
evalList = traverse . eval

env :: Data.Map.Map String Int
env =
  Data.Map.fromList
    [ ("x", 5),
      ("y", 6)
    ]

exps =
  [ sampleExprWithAddSubMulDivAndVars,
    evenBiggerSampleExprWithAddSubMulDivAndVars
  ]

-- >>> evalList' env exps
-- [Right 0,Left "Variable z is missing."]
-- >>> evalList env exps
-- Left "Variable z is missing."

simpleProgramWithLet =
  Let
    "x"
    (Val 5)
    ( Add
        (Var "x")
        ( Let
            "y"
            (Add (Val 1) (Val 2))
            ( Add (Var "x") (Var "y")
            )
        )
    )

-- >>> eval Data.Map.empty simpleProgramWithLet
-- Right 13
