(* Kolejki dwustronne *)
(* Implementacja za pomocą list i produktów kartezjańskich. *)

open List;;
#use "fifo.ml";;

module type BIQUEUE = 
  sig
    include QUEUE
    val back : 'a queue -> 'a
    val insert_front : 'a queue -> 'a -> 'a queue
    val remove_back : 'a queue -> 'a queue
    val size : 'a queue -> int
  end;;

module FifoLifo : BIQUEUE =
  struct
    exception EmptyQueue

    (* Typ = przód, tył, rozmiar. *)
    type 'a queue = 'a list * 'a list * int

    (* Pusta kolejka *)
    let empty : 'a queue = ([], [], 0)
	
    (* Rozmiar *)
    let size ((_, _, n): 'a queue) = n

    (* Czy jest pusta? *)
    let is_empty q = size q = 0
 
    (* Wyrównanie listy, jeśli jedna z połówek jest pusta, a druga 
       zawiera więcej niż jeden element. *)
    let balance (q: 'a queue) = 
      (* Podział listy na połowy h1 i h2, length h1 <= lenght h2 *)
      let halve l = 
      	let rec head l n a = 
      	  if n = 0 then 
                  (rev a, l)
      	  else
                  match l with 
                    []   -> failwith "Zbyt krótka lista" |
                    h::t -> head t (n-1) (h::a)
      	in 
        	head l (length l / 2) []
      in 
        match q with
        	([], [], _)    -> q |
        	([_], [], _)   -> q |
        	([], [_], _)   -> q |
        	(front, [], n) -> 
            let (h1, h2) = halve front
            in (h1, rev h2, n) |
        	([], back, n)  -> 
            let (h1, h2) = halve back
            in (rev h2, h1, n) |
    	    _           -> q 
        		
    (* Wstawienie na koniec *)
    let insert ((front, back, n): 'a queue) x =
      ((balance (front, x::back, n+1)): 'a queue)

    (* Wstawienie na początek *)
    let insert_front ((front, back, n): 'a queue) x =
      ((balance (x::front, back, n+1)): 'a queue)

    (* Początek *)
    let front (q: 'a queue) = 
      match q with
      	([], [x], _) -> x |
      	([], _, _)   -> raise EmptyQueue |
      	(x::_, _, _) -> x

    (* Koniec *)
    let back (q: 'a queue) = 
      match q with
      	([x], [], _) -> x |
      	(_, [], _)   -> raise EmptyQueue |
      	(_, x::_, _) -> x
  
    (* Wyjęcie początku *)
    let remove (q: 'a queue) =
      match q with
      	(_::front, back, n) -> balance (front, back, n-1) |
        ([], [_], _) -> empty |
      	([], _, _) -> raise EmptyQueue

    (* Wyjęcie końca *)
    let remove_back (q: 'a queue) =
      match q with
      	(front, _::back, n) -> balance (front, back, n-1) |
        ([_], [], _) -> empty |
      	(_, [], _) -> raise EmptyQueue
    
  end;;



(* Testy *)
open FifoLifo;;

assert (try ignore (front empty); false with EmptyQueue -> true);;
assert (try ignore (back empty); false with EmptyQueue -> true);;
assert (try ignore (remove empty); false with EmptyQueue -> true);;
assert (try ignore (remove_back empty); false with EmptyQueue -> true);;

let q1 = insert (insert (insert empty 1) 2) 3;;
assert (front q1 = 1);;
assert (back q1 = 3);;
let q2 = remove_back (remove q1);;
assert (front q2 = back q2);;
let q3 = remove_back q2;;
assert (size q3 = 0);;