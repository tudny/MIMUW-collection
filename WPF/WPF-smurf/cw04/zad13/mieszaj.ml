(* z zadania 6 *)

let rec shuffle l1 l2 =
    match (l1, l2) with
    | (_, []) -> l1
    | ([], _) -> l2
    | (h1 :: t1, h2 :: t2) -> h1 :: h2 :: (shuffle t1 t2);;

let devide lista =
  let rec aux (h :: t) ile acc =
    if ile = 0 then (List.rev acc, h :: t)
    else aux t (ile - 1) (h :: acc) in
  let cnt = (1 + List.length lista) / 2 in
  match lista with
  | []     -> ([], [])
  | [x]    -> ([x], [])
  | _ :: _ -> aux lista cnt [];;

(* devide [];; *)
devide [1];;
devide [1; 2];;
devide [1; 2; 3];;
devide [1; 2; 3; 4];;
devide [1; 2; 3; 4; 5];;
devide [1; 2; 3; 4; 5; 6];;


let mieszaj lista =
  let (first, second) = devide lista in
  shuffle first (List.rev second);;

mieszaj [1; 2; 3; 4; 5];;
mieszaj [1; 2];;
mieszaj [1; 2; 3];;
mieszaj [1; 2; 3; 4];;
