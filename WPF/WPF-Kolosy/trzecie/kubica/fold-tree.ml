(* Drzewo binarne. *)
type 'a tree = 
  Node of 'a * 'a tree * 'a tree | 
  Leaf;;

(* Fold w górê drzewa. *)
let rec fold_tree f a t = 
  match t with
    Leaf -> a |
    Node (x, l, r) -> f x (fold_tree f a l) (fold_tree f a r);;

(* Zadanie o liczbie widocznych elementów. *)
let widoczne t = 
  let merge x l r = 
    fun k -> if x < k then l k + r k else l x + r x + 1
  in
    fold_tree merge (fun _ -> 0) t (-1);;

assert (widoczne (Node (6,
			Node (4, Leaf, Leaf), 
			Node (8, 
			      Node (10, Leaf, Leaf), 
			      Node (5, Leaf, Leaf))
		       )
		 )=3);;

(* Zadanie o odchyleniach i szeroko¶ci drzewa. *)
let szerokosc t = 
  let f _ (minl, maxl) (minr, maxr) = (min (minl-1) (minr+1), max (maxl-1) (maxr+1))
  in 
    let (x, y) = fold_tree f (1,-1) t
    in y-x;;

assert (szerokosc (Node (42, Leaf, Leaf)) = 0);;
assert (szerokosc (Node (6,
			Node (4, Leaf, Leaf), 
			Node (8, 
			      Node (10, 
				    Leaf, 
				    Node (5, 
					  Leaf, 
					   Node (9, Leaf, Leaf))), 
			      Leaf)
		       )
		 )=3);;


open List;;

(* Drzewa dowolnego stopnia. *)
type 'a tree =  Node of 'a * 'a tree list;;

(* fold_tree *)
let rec fold_tree f (Node (x, l)) = 
  f x (map (fold_tree f) l);;


(* Wysoko¶æ drzewa *)
let height t = 
  let maxl l = fold_left max 0 l
  in fold_tree (fun _ l -> maxl l + 1) t;;

assert (height (Node ((),[])) = 1);;
assert (height (Node (1,[Node (2,[Node (4,[])]); Node (3,[])])) = 3);;


(* map_tree *)
let map_tree f t = 
  fold_tree (fun x l -> Node (f x, l)) t;;

assert (map_tree (fun x -> x+1) (Node (1,[Node (2,[Node (4,[])]); Node (3,[])])) = Node (2, [Node (3, [Node (5, [])]); Node (4, [])]));;


(* Zadanie o liczbie widocznych elementów w drzewie. *)
let widoczne t = 
  let merge x l k = 
    if x < k then 
      fold_left (fun a h -> a + h k) 0 l
    else 
      fold_left (fun a h -> a + h x) 0 l + 1
  in
    fold_tree merge t (-1);;

assert (widoczne (Node (6, [
			Node (4, []);
			Node (8, [
			      Node (10, []);
			      Node (5, []);
			      Node (7, [
				    Node (12, [])
				  ])
			    ])
		      ])
		 )=4);;

