
let fit c lista = 
  let rec aux l1 l2 acc =
    match l1 with
    | [] -> acc
    | [x] -> acc
    | (h1 :: t1) -> 
      match l2 with 
      | [] -> aux t1 (List.tl t1) acc
      | (h2 :: t2) ->
        let nowy_acc = min acc (abs (c + h1 + h2)) in
        aux l1 t2 nowy_acc
  in aux lista (List.tl lista) max_int;;

fit 1 [-42; -12; -8; -1; -1; 5; 15; 60];;
fit 2 [1; 2];;
fit 3 [-1; 2; 3];;
fit 4 [-100; 0; 200];;
fit 42 [-28; -25; -15; -1; 4; 8; 15; 60];;
