
let rec shuffle l1 l2 =
    match (l1, l2) with
    | (_, []) -> l1
    | ([], _) -> l2
    | (h1 :: t1, h2 :: t2) -> h1 :: h2 :: (shuffle t1 t2);;

shuffle [3; 2; 8; 1; 9; 3; 6] [5; 7; 0];;
shuffle [5; 7; 0] [3; 2; 8; 1; 9; 3; 6];;
shuffle [1; 2; 3] [4; 5; 6];;
