
open List;;

type 'a tree =  Node of 'a * 'a tree list;;

let rec fold_tree f (Node (x, l)) =
  f x (map (fold_tree f) l);;


let preoder (Node(v, li) as tr) =
  let process v li acc =
    List.fold_left
      (fun _acc h -> h _acc)
      (v :: acc)
      li
  in List.rev (fold_tree process tr [])
;;


let t = Node (1, [Node(2, [Node(4, []); Node(5, [])]); Node(3, [])]);;
preoder t;;


let postorder (Node(v, li) as tr) =
  let process v li acc =
    v :: List.fold_left
      (fun _acc h -> h _acc)
      (acc)
      li
  in List.rev (fold_tree process tr [])
;;

let t = Node (1, [Node(2, [Node(4, []); Node(5, [])]); Node(3, [])]);;
postorder t;;
