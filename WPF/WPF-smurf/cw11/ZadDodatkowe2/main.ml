
type 'a tree = Node of 'a * 'a tree list;;

let rec fold_tree f (Node (x, l)) =
  f x (List.map (fold_tree f) l);;


(* Rozwiązanie O(n log n) *)

let parse li =
  let rec loop l akt_wys acc =
    print_int akt_wys; print_string "\n";
    match l with
    | [] -> acc
    | h :: t ->
      let (a_h :: a_t) = acc in
        if h = akt_wys then loop t akt_wys (a_h + 1 :: a_t)
        else loop l (akt_wys + 1) (0 :: acc)
  in loop li 0 [0];;

let liscie t =
  let process v li (h, acc) =
    match li with
    | [] -> (h :: acc)
    | _  ->
      List.fold_left
        (fun acc f -> f (h + 1, acc))
        acc li
  in fold_tree process t (0, []) |>
                        List.rev |>
                List.sort compare|>
                           parse |> List.rev
;;


let t = Node(1, [Node(2, [Node(3, []); Node(4, [])]); Node(5, [Node(6, [Node(7, [Node(9, [])]); Node(8, [])])])]);;
liscie t;;


