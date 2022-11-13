
(* Reprezetnacja grafu : tablica list sąsiedztwa *)

type graph = int list array;;

let path g = 
  let n = Array.length g in 
  let dp = Array.make n (-1) in
  Array.iteri ( fun i edges -> 
    dp.(i) <- 1 +
      List.fold_left ( fun acc e -> 
        max acc dp.(e)
      ) (-1) edges
    ) g;
  Array.fold_left (max) (-1) dp
;;

let g = 
  [|
    [1; 2; 3];
    [0; 3; 5];
    [0; 3];
    [0; 1; 2; 4; 5];
    [3; 5];
    [1; 3; 4; 7; 8];
    [5];
    [5; 8];
    [5; 7]
  |]
;;

Printf.printf "path: %d" (path g);;
