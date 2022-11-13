open List;;
#use "listy.ml";;
#use "pierwiastki.ml";;



(* Selection sort *)

(* Podzia³ listy na element maksymalny i resztê. *)
let select_max (h::t) =
  fold_left 
    (fun (m, r) x -> if x > m then (x, m::r) else (m, x::r)) 
    (h, []) t;;

(* Sortowanie *)
let selection_sort l =
  let przenies (s, l) = 
    let (m, r) = select_max l
    in (m::s, r)
  in 
    iteruj 
      ([], l) 
      przenies 
      (fun (_, l) -> l = []) 
      (fun (s, _) -> s);;


(* Insertion sort *)
let wstaw l x = 
  (filter (fun y -> y <= x) l) @ (x :: (filter (fun y -> y > x) l));;

let insertion_sort l =  fold_left wstaw [] l;;


(* Insertion sort - koszt O(n + liczba inwersji). *)
let wstaw x l = 
  let rec wst a x l = 
    match l with 
      []   -> rev(x::a) | 
      h::t -> 
	if x > h then wst (h::a) x t 
	else rev (x::a) @ l
  in wst [] x l;;

let insertion_sort l =  fold_right wstaw l [];;


(* Merge sort *)

let rec merge_sort l =
  let split l = 
    fold_left (fun (l1, l2) x -> (l2, x::l1)) ([], []) l
  and merge l1 l2 = 
    let rec mrg a l1 l2 = 
      match l1 with 
        [] -> (rev a) @ l2 |
        (h1::t1) -> 
          match l2 with 
            [] -> (rev a) @ l1 |
            (h2::t2) -> 
              if h1 > h2 then 
                mrg (h2::a) l1 t2
              else
                mrg (h1::a) t1 l2
    in 
      mrg [] l1 l2
  in 
    match l with
      [] -> []|
      [x] -> [x]|
      _ -> 
        let 
          (l1, l2) = split l
        in 
          merge (merge_sort l1) (merge_sort l2);;

(* Quick sort *)
let rec quick_sort l = 
  let split l x = (filter (fun y -> y < x) l, 
                   filter (fun y -> y = x) l, 
                   filter (fun y -> y > x) l) 
  in 
    if length l < 2 then l else
      let  x = nth l (Random.int (length l))
      in
        let (ll, le, lg) = split l x
        in
          (quick_sort ll) @ le @ (quick_sort lg);;


(* Kolejka priorytetowa *)
module type PRI_QUEUE = sig
  type 'a pri_queue
  val empty_queue : 'a pri_queue
  val is_empty : 'a pri_queue -> bool
  val put : 'a pri_queue -> 'a -> 'a pri_queue
  val getmax : 'a pri_queue -> 'a 
  val removemax : 'a pri_queue -> 'a pri_queue
  exception Empty_Queue
end;;

(* Kolejka priorytetowa jako lista nieuporz±dkowana *)
module Unordered_Pri_Queue : PRI_QUEUE = struct
  exception Empty_Queue
  type 'a pri_queue = 'a list
  let empty_queue = []
  let is_empty q = q = []
  let put q x = x::q
  let getmax q =
    if q = [] then raise Empty_Queue
    else fst (select_max q)
  let removemax q =
    if q = [] then raise Empty_Queue
    else snd (select_max q)
end;;


(* Kolejka priorytetowa jako lista uporz±dkowana nierosn±co *)
module Ordered_Pri_Queue : PRI_QUEUE = struct
  exception Empty_Queue
  type 'a pri_queue = 'a list
  let empty_queue = []
  let is_empty q = q = []
  let put q x = 
    (filter (fun y -> y > x) q) @ 
    (x :: (filter (fun y -> y <= x) q))
  let getmax q =
    if q = [] then raise Empty_Queue
    else hd q
  let removemax q =
    if q = [] then raise Empty_Queue
    else tl q
end;;

(* Kolejka priorytetowa jako heap *)
module Heap_Pri_Queue : PRI_QUEUE = struct
  exception Empty_Queue

  type 'a pri_queue =
    Node of 'a * 'a pri_queue * 'a pri_queue * int |
    Leaf 

  let empty_queue = Leaf

  let is_empty q = q = Leaf

  (* Rozmiar *)
  let size q = 
    match q with
      Leaf -> 0 |
      Node (_, _, _, n) -> n

  (* Pomocniczy selektor. *)
  let getmax h = 
    match h with 
      Leaf -> raise Empty_Queue |
      Node (r, _, _, _) ->  r

  (* Pomocniczy modyfikator. *)
  let set_root h r = 
    match h with 
      Leaf -> Node (r, Leaf, Leaf, 1) |
      Node (_, l, p, n) -> Node (r, l, p, n) 

  let rec put h x =
    match h with
      Leaf -> Node (x, Leaf, Leaf, 1) |
      Node (r, l, p, n) ->
        if size l <= size p then 
	  Node((max x r), (put l (min x r)), p, (n+1))
        else 
	  Node((max x r), l, (put p (min x r)), (n+1))
  
  let rec removemax h =
      match h with
	Leaf -> raise Empty_Queue |
	Node (_, Leaf, Leaf, _) -> Leaf |
	Node (_, l, Leaf, _) -> l |
	Node (_, Leaf, p, _) -> p |
	Node (_, (Node (rl, _, _, _) as l),
                 (Node (rp, _, _, _) as p), n) ->
	  if rl >= rp then
            Node (rl, removemax l, p, n - 1)
          else
            Node (rp, l, removemax p, n - 1)
end;;


(* Sortowanie w oparciu o kolejkê priorytetow± *)
module type HEAP_SORT = functor (Queue : PRI_QUEUE) -> sig
  val heap_sort : 'a list -> 'a list
end;;

module Heap_Sort : HEAP_SORT = functor (Queue : PRI_QUEUE) -> struct
  open Queue
  let heap_sort l = 
    let wloz = fold_left put empty_queue l
    and wyjmij (l, q) = ((getmax q)::l, removemax q)
    in
      iteruj
        ([], wloz)
        wyjmij
        (fun (_, q) -> is_empty q)
        (fun (l, _) -> l)
end;;

module Sort = Heap_Sort(Unordered_Pri_Queue);;
Sort.heap_sort [1;4;2;8;3;12;5;9;6;11;7];;

module Sort = Heap_Sort(Ordered_Pri_Queue);;
Sort.heap_sort [1;4;2;8;3;12;5;9;6;11;7];;

module Sort = Heap_Sort(Heap_Pri_Queue);;
Sort.heap_sort [1;4;2;8;3;12;5;9;6;11;7];;
