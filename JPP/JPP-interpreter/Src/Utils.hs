
module Src.Utils where

left :: (b -> c) -> Either b d -> Either c d
left m (Left x) = Left (m x)
left _ (Right r) = Right r
