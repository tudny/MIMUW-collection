
let rec rem_dup li =
    match li with
    | []     -> []
    | [x]    -> [x]
    | h1 :: h2 :: t ->
        (if h1 = h2 then [] else [h1]) @ (rem_dup (h2 :: t));;

rem_dup [1; 2; 3; 3; 4; 4];;
rem_dup [1; 2; 3; 4];;
rem_dup [1; 2; 2; 3; 4];;
rem_dup [1; 1; 2; 3; 4; 1];;
rem_dup [];;
rem_dup [1];;


let count_distinct l1 l2 =
    let a = rem_dup l1
    and b = rem_dup l2 in
        let rec aux n m cnt =
            match (n, m) with
            | (_, []) -> cnt + List.length n
            | ([], _) -> cnt + List.length m
            | (h1 :: t1, h2 :: t2) ->
                if h1 = h2 then aux t1 t2 (cnt + 1)
                else if h1 < h2 then aux t1 m (cnt + 1)
                else (* if h1 > h2 *) aux n t2 (cnt + 1)
        in aux a b 0;;


count_distinct [1; 2; 3; 4; 4; 5; 6] [3; 4; 5; 10; 11];;
count_distinct [3; 4; 5; 10; 11] [1; 2; 3; 4; 4; 5; 6];;
count_distinct [1; 2; 3; 4; 4; 5; 6] [1; 2; 3; 4; 4; 5; 6];;
count_distinct [] [];;
count_distinct [1; 2; 3] [2; 3; 4];;
