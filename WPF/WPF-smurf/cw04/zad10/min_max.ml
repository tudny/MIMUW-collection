
let min_max lis =
  let update_kandydaci kand ele =
    let (mini, maks) = kand in
    if ele < mini then (ele, maks)
    else if ele > maks then (mini, ele)
    else kand
  in
  let rec aux lista kandydaci =
    match lista with
    | []     -> kandydaci
    | [x]    -> update_kandydaci kandydaci x
    | h1 :: h2 :: t ->
            let (mini, maks) = kandydaci in
            if h1 > h2 then aux t (min mini h2, max maks h1)
            else aux t (min mini h1, max maks h2)
  in match lis with
  | [] -> (max_int, min_int)
  | (h :: t) -> aux t (h, h);;

min_max [1; 2; 3; 4; 5; 6; 7; 8];;
min_max [10; 5; 15; 0; 20; (-5); 25; (-10)];;
