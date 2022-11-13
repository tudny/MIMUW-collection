#use "listy.ml";;
#use "fifo.ml";;


(* Drzewo kodowania Huffmana *)
type 'a huff_tree = 
  Letter of ('a * float) |
  Node of ('a huff_tree * 'a huff_tree * float);;

let frequency (t : 'a huff_tree) = 
  match t with
    Letter (_, f) -> f |
    Node (_, _, f)   -> f;;

let merge t1 t2 = Node (t1, t2, frequency t1 +. frequency t2);;


(* Kolekcja drzew Huffmana *)
type 'a huff_col = 'a huff_tree fifo * 'a huff_tree fifo;; 


(* Tworzy kolekcje lisci Huffmana na podstawie listy czestosci *)
(* wystepowania uporzadkowanej niemalejaco.                    *)
let make_huff_col l =
  (fold (fun q x -> put q (Letter x)) empty_queue l, empty_queue);;

let huff_size ((q1, q2): 'a huff_col) = size q1 + size q2;;

let huff_put ((q1, q2): 'a huff_col) t = ((q1, put q2 t): 'a huff_col);;
  
let huff_first ((q1, q2): 'a huff_col) = 
  if q1 = empty_queue then first q2 else 
  if q2 = empty_queue then first q1 else
    let f1 = first q1
    and f2 = first q2
    in 
      if frequency f1 <= frequency f2 then f1 else f2;;

let huff_remove ((q1, q2): 'a huff_col) =
  if q1 = empty_queue then (q1, remove q2) else 
  if q2 = empty_queue then (remove q1, q2) else
    let f1 = first q1
    and f2 = first q2
    in
      if frequency f1 <= frequency f2 then 
        (remove q1, q2) 
      else 
        ((q1, remove q2): 'a huff_col);;

(* Algorytm Huffmana *)
let huffman l = 
  let rec process col = 
    if huff_size col = 1 then 
      huff_first col
    else 
      let t1 = huff_first col
      and col1 = huff_remove col
      in 
        let t2 = huff_first col1
        and col2 = huff_remove col1
        in 
          process (huff_put col2 (merge t1 t2))
  in 
    process (make_huff_col l);;
  
(* Przyk³adowe dane *)
huffman ['C',1.; 'D',1.; 'E',1.; 'F',1.; 'G',1.; 'H',1.; 'B',3.; 'A',8.];;
(* = 
Node
 (Letter ('A', 8),
  Node
   (Node
     (Node (Letter ('C', 1), Letter ('D', 1), 2),
      Node (Letter ('E', 1), Letter ('F', 1), 2), 4),
    Node (Node (Letter ('G', 1), Letter ('H', 1), 2), Letter ('B', 3), 5), 9),
  17)
*)
  
  
  
  
