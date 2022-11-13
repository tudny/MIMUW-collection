(* Kolejki fifo *)
(* Implemetacja za pomoca list i rekordów *)

open List;;


module type QUEUE = 
  sig
    exception EmptyQueue
    type 'a queue
    val empty : 'a queue
    val is_empty : 'a queue -> bool
    val insert : 'a queue -> 'a -> 'a queue
    val front : 'a queue -> 'a
    val remove : 'a queue -> 'a queue
  end;;

module Lifo : QUEUE = 
  struct
    exception EmptyQueue
    type 'a queue = 'a list
    let empty = []
    let is_empty q = q = []
    let insert q x = x::q
    let front q = 
      if q = [] then raise EmptyQueue 
      else hd q
    let remove q =
      if q = [] then raise EmptyQueue 
      else tl q
  end;;
	  
module Fifo : QUEUE = 
  struct
    exception EmptyQueue

    (* Kolejka to trójka: przód, ty³, rozmiar. *)
    (* Je¿eli przód jest pusty, to kolejka jest pusta. *)
    type 'a queue = {front: 'a list; back: 'a list; size: int}

    let empty = {front=[]; back=[]; size=0}

    let size q = q.size

    let is_empty q = size q = 0

    let balance q = 
      match q with
	{front=[]; back=[]}        -> q |
	{front=[]; back=b; size=s} -> {front=rev b; back=[]; size=s} |
	_                          -> q

    let insert {front=f; back=b; size=n} x =
      balance {front=f; back=x::b; size=n+1}

    let front q = 
      match q with
	{front=[]}   -> raise EmptyQueue |
	{front=x::_} -> x
    
    let remove q =
      match q with
	{front=_::f} -> balance {front=f; back=q.back; size=q.size-1} |
	_            -> raise EmptyQueue

  end;;

