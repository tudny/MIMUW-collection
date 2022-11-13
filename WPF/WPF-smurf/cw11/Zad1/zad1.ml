
open List;;

type 'a tree =  Node of 'a * 'a tree list;;

let rec fold_tree f (Node (x, l)) =
  f x (map (fold_tree f) l);;

let zlicz t li =
  fold_left (+) 0 li + 1;;

let ile_wezlow t =
  fold_tree zlicz t;;

let t1 = Node(10, []);;
ile_wezlow t1;;

let t2 = Node(10, [Node(20, [])]);;
ile_wezlow t2;;
