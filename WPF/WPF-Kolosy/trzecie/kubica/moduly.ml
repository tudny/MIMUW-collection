#use "listy.ml";;

(* Przyk³ady prostych struktur. *)

module Modulik = 
  struct
    type typik = int list
    let lista = [2; 3; 7]
    let rec prod l = 
      if l = [] then 1 else hd l * prod (tl l) 
  end;;

Modulik.prod Modulik.lista;;
open Modulik;;
prod lista;;

module Pusty = 
  struct
  end;;

module M = 
  struct
    module A = 
      struct
	let a = 27
      end
    module B = 
      struct 
	let b = 15
      end
  end;;

M.A.a + M.B.b;;


(* Przyk³ady sygnatur *)

module type S = 
  sig
    type abstrakcyjny
    type konkretny = int * float
    val x : abstrakcyjny * konkretny
    module type Pusty = sig end
    exception Wyjatek of abstrakcyjny
  end;;

module type FIFO = 
  sig
    exception EmptyQueue
    type 'a queue
    val empty : 'a queue
    val insert : 'a queue -> 'a -> 'a queue
    val front : 'a queue -> 'a
    val remove : 'a queue -> 'a queue
  end;;


(* Przyk³ady struktur *)

module Empty : sig end = struct end;;


module Fifo_implementation = 
  struct
    exception EmptyQueue
    type 'a queue = 'a list
    let empty = []
    let insert q x = q @ [x]
    let front q = 
      match q with
	[] -> raise EmptyQueue |
	h::_ -> h
    let remove q = 
      match q with
	[] -> raise EmptyQueue |
	_::t -> t
  end;;

module Fifo : FIFO = Fifo_implementation;;
module Fifo = (Fifo_implementation : FIFO);;

module Fifo : 
    sig
      exception EmptyQueue
      type 'a queue
      val empty : 'a queue
      val insert : 'a queue -> 'a -> 'a queue
      val front : 'a queue -> 'a
      val remove : 'a queue -> 'a queue
    end = 
  struct
    exception EmptyQueue
    type 'a queue = 'a list
    let empty = []
    let insert q x = q @ [x]
    let front q = 
      match q with
	[] -> raise EmptyQueue |
	h::_ -> h
    let remove q = 
      match q with
	[] -> raise EmptyQueue |
	_::t -> t
  end;;

module type QUEUE = 
  sig 
    include FIFO 
    val back : 'a queue -> 'a
    val insert_front : 'a queue -> 'a -> 'a queue
    val remove_back : 'a queue -> 'a queue
  end;;
	
module type LIST =  FIFO with type 'a queue = 'a list;;

module Q = (Fifo_implementation : LIST);;


(* Liczby wymierne, bez funktorów *)
module type ULAMKI = 
  sig 
    type t
    val ulamek : int -> int -> t
    val licznik : t -> int
    val mianownik : t -> int
  end;;

module type RAT = 
  sig
    include ULAMKI 
    val plus : t -> t -> t
    val minus : t -> t -> t
    val razy : t -> t -> t
    val podziel : t -> t -> t
    val rowne : t -> t -> bool
  end;;

module Ulamki : ULAMKI = 
  struct
    type t = int * int
    let ulamek l m = (l, m)
    let licznik (l, _) = l
    let mianownik (_, m) = m
  end;;

module Rat : RAT =
    struct
      include Ulamki 
      let plus x y = 
	ulamek 
	  (licznik x * mianownik y + licznik y * mianownik x) 
	  (mianownik x * mianownik y)
      let minus x y = 
	ulamek 
	  (licznik x * mianownik y - licznik y * mianownik x) 
	  (mianownik x * mianownik y)
      let razy x y = 
	ulamek
	  (licznik x * licznik y)
	  (mianownik x * mianownik y)
      let podziel x y = 
	ulamek
	  (licznik x * mianownik y)
	  (mianownik x * licznik y)
      let rowne x y = 
	(licznik x * mianownik y) = (licznik y * mianownik x)
    end;;


