
open List;;

type 'a tree =  Node of 'a * 'a tree list;;

let rec fold_tree f (Node (x, l)) =
  f x (map (fold_tree f) l);;

let int_of_bool b = if b then 1 else 0

(* acc (ile_wierz, ile_dobrze) *)
let is_regular (Node (v, li) as tr) =
  if li = [] then true else
  let cnt = length li in
  let regul t li =
    match li with
    | [] -> (1, 1)
    | _ -> List.fold_left
            (fun (ile, tru) (ile_w, ile_true) -> (ile + ile_w, tru + ile_true))
            (1, int_of_bool (List.length li + 1 = cnt)) li
  in let (ile_w_sumie, ile_dobrych) = fold_tree regul tr in
  ile_w_sumie - 1 = ile_dobrych
;;


let t1 = Node(1, [Node(2, [])]);;
is_regular t1;;

let t1 = Node(1, []);;
is_regular t1;;

let t2 = Node(1, [Node(2, []); Node(3, [])]);;
is_regular t2;;

let t2 = Node(1, [Node(2, [Node(4, [])]); Node(3, [Node(5, [])])]);;
is_regular t2;;

let t2 = Node(1, [Node(2, [Node(4, []); Node(6, [])]); Node(3, [Node(5, [])])]);;
is_regular t2;;


(* Syyyyf *)
(*

let is_regualar (Node(v, li) as tr) =
  let cnt = length li + 1 in
  let rec add n acc = if n < 1 then acc else add (n - 1) ((Node(0, [])) :: acc) in
  let create_empty = add (cnt - 2) [tr] in
  let super_tr = Node(0, create_empty) in
  let regualar t li =
    match li with
    | [] -> true
    | _ :: _ -> fold_left (&&) (length li + 1 = cnt) li
  in
  fold_tree (regualar) super_tr
;;

let t1 = Node(1, [Node(2, [])]);;

is_regualar t1;; *)






