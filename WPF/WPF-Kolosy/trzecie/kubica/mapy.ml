(* M A P Y *)


(* Sygnatura mapy. *)
module type MAP = 
  sig
    (* Mapa z warto¶ci typu 'a w 'b. *)
    type ('a,'b) map 

    (* Wyj±tek podnoszony, gdy badamy warto¶æ spoza dziedziny. *)
    exception Undefined

    (* Pusta mapa. *)
    val empty : ('a,'b) map
	
    (* Predykat charakterystyczny dziedziny mapy. *)
    val dom : ('a,'b) map -> 'a -> bool

    (* Zastosowanie mapy. *)
    val apply : ('a,'b) map -> 'a -> 'b

    (* Dodanie warto¶ci do mapy. *)
    val update : ('a,'b) map -> 'a -> 'b -> ('a,'b) map
  end;;


  
(* Bardzo prosta mimplementacja proceduralna map. *)
module ProcMap : MAP = 
  struct
    (* Mapa z warto¶ci typu 'a w 'b. *)
    type ('a, 'b) map = 'a -> 'b

    (* Wyj±tek podnoszony, gdy badamy warto¶æ spoza dziedziny. *)
    exception Undefined

    (* Pusta mapa *)
    let empty = function _ -> raise Undefined

    (* Predykat charakterystyczny dziedziny mapy. *)
    let dom m x = 
      try 
	let _ = m x 
	in true
      with Undefined -> false

    (* Zbadanie warto¶ci w punkcie *)
    let apply m x = m x

    (* Dodanie warto¶ci do mapy. *)
    let update f a v = 
      function x ->
	if a = x then v else f x
  end;;

(*********************************************************************)

(* Mapy jako drzewa BST *)
module BstMap : MAP =
  struct
    (* Mapa z warto¶ci typu 'a w 'b. *)
    type ('a, 'b) map = Empty | Node of ('a * 'b * ('a, 'b) map * ('a, 'b) map)

    (* Wyj±tek podnoszony, gdy badamy warto¶æ spoza dziedziny. *)
    exception Undefined

    (* Pusta mapa to puste drzewo *)   
    let empty = Empty

    (* Znajd¼ poddrzewo o zadanym kluczu              *)
    (* (je¶li nie ma, to wynikiem jest puste drzewo). *)
    let rec find tree key =
      match tree with
	Empty -> Empty |
	Node (k, _, l, r) ->
	  if key = k then tree
	  else if key < k then find l key
	  else find r key

    (* Zastosowanie mapy. *)
    let apply m k =
      match find m k with
	Empty -> raise Undefined |
	Node (_, v, _, _) -> v
    
    (* Sprawdzenie, czy punkt nale¿y do dziedziny *)
    let dom m x =
      try
	let _ = apply m x
	in true
      with
	Undefined -> false

    (* Zmiana warto¶ci mapy w punkcie *)
    let rec update tree key value =
      match tree with
	Empty -> Node (key, value, Empty, Empty) |
	Node (k, v, l, r) ->
	  if key = k then
            Node (key, value, l, r)
	  else if key < k then
            Node (k, v, update l key value, r)
	  else
            Node (k, v, l, update r key value)
  end;;
