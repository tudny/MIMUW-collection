let rec print_list = function
    [] -> print_endline ""
| e::l -> print_int e ; print_string ", " ; print_list l


let rec make_k_empty k =
  if k = 0 then []
  else [] :: (make_k_empty (k - 1));;

let podziel lista k =
  let rec aux left left_to_add acc =
    print_list left;
    match left_to_add with
    | [] -> aux left (List.rev acc) []
    | (head :: tail) ->
      match left with
      | [] -> acc
      | (h :: t) ->
        aux t tail ((h :: head) :: acc)
  in aux lista (make_k_empty k) [];;

podziel [1; 2; 3; 4; 5; 6; 7] 3;;
