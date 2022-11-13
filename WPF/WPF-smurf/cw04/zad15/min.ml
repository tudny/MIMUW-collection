
let malo lista =
  let rec aux l1 l2 acc =
    match l1 with
    | [] -> acc
    | [x] -> acc
    | (h1 :: t1) ->
      match l2 with
      | [] -> aux t1 (List.tl t1) acc
      | (h2 :: t2) ->
        let nowy_acc = min acc (abs (h1 + h2)) in
        aux l1 t2 nowy_acc
  in aux lista (List.tl lista) max_int;;

malo [-42; -12; -8; -1; -1; 5; 15; 60];;
malo [1; 2];;
malo [-1; 2; 3];;
malo [-100; 0; 200];;
