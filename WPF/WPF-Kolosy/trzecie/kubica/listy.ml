(* Rekurencyjne implementacje niektórych operacji na listach. *)
open List;;

(* n-ty element *)
let rec nth l n = 
  if n = 0 then hd l 
  else nth (tl l) (n - 1);;

(* Odwrócenie listy *)
let rev l =
  let rec pom l w = 
    if l = [] then w
    else pom (tl l) ((hd l) :: w)
  in pom l [];;

(* Sklejanie list. *)
let rec append l1 l2 =
    if l1 = [] then l2
    else (hd l1) :: (append (tl l1) l2);;

    
(* Przydatne procedury na listach *)  
let rec fold_right f l a = 
  match l with
    []   -> a |
    h::t -> f h (fold_right f t a);;

let rec fold_right2 f l1 l2 a = 
  match (l1, l2) with 
    ([], [])         -> a |
    (h1::t1, h2::t2) -> f h1 h2 (fold_right2 f t1 t2 a) |
    _                -> failwith "Listy ró¿nej d³ugo¶ci";;

let rec fold_left f a l = 
  match l with 
    []   -> a |
    h::t -> fold_left f (f a h) t;;

let rec fold_left2 f a l1 l2 = 
  match (l1, l2) with
    ([], [])         -> a |
    (h1::t1, h2::t2) -> fold_left2 f (f a h1 h2) t1 t2 |
    _                -> failwith "Listy ró¿nej d³ugo¶ci";;

let iloczyn_skalarny l1 l2 = 
  fold_left2 (fun a x y -> a + x * y) 0 l1 l2;;

let map f l = fold_right (fun h t -> (f h)::t) l [];;

let map2 f l1 l2 = fold_right2 (fun h1 h2 t -> (f h1 h2)::t) l1 l2 [];;

let suma_wektorow l1 l2 = map2 (+) l1 l2;;

let filter p l = 
  fold_right (fun h t -> if p h then h::t else t) l [];;

let rev l = fold_left (fun v x -> x::v) [] l;;

let append l1 l2 = fold_right (fun x v -> x::v) l1 l2;;

let flatten l = fold_right (@) l [];;


(* Jeden fold za pomoc± drugiego. *)
let fold_right f l a = 
  let rev = fold_left (fun a h -> h::a) [] 
  in fold_left (fun x h -> f h x) a (rev l);;

let fold_left f a l = 
  fold_right (fun h p -> function x -> p (f x h)) l (fun x -> x) a;;

(* Sito Eratostenesa *)
let rec ints a b = if a > b then [] else a :: ints (a+1) b;;

let rec sito l = 
  match l with
    [] -> []|
    h::t -> h::sito (filter (function x -> x mod h <> 0) t);;

let eratostenes n = sito (ints 2 n);;

(* Przekszta³ca parê list w listê par -- odpowiednik produktu. *)
let pairs l1 l2 = 
  flatten (map (fun i-> map (fun j -> (i,j)) l1) l2);;

(* Lista g³ów list sk³adowych. *)
let heads l = 
  flatten (map (function [] -> [] | h::_ -> [h]) l) ;;


(* Czy istnieje element spe³niaj±cy predykat? *)
exception Exists;;

let exists p l = 
  try fold_left (fun a x -> if p x then raise Exists else a) false l
  with Exists -> true;;


