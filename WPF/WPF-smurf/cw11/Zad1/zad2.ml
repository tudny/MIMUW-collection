
let non f x = not (f x);;

let exists p li =
  List.fold_left (fun acc el -> acc || p el) false li;;

let forall p li =
  not (exists (non p) li);;

open List;;

type 'a tree =  Node of 'a * 'a tree list;;

exception Exit;;

let rec fold_tree f (Node (x, l)) =
  f x (map (fold_tree f) l)
;;

let rownowaga v li =
  let wys =
    try List.hd li
    with | _ -> 0
  in
  if forall (fun h -> h = wys) li then wys + 1
  else raise Exit
;;

let czy_rownowazne t =
  try
    let _ = fold_tree rownowaga t in true
  with
  | Exit -> false


let t1 = Node(1, [Node(2, []); Node(3, [Node(4, [])])]);;

czy_rownowazne t1;;


