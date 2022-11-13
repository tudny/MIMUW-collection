open List;;
(* Czê¶ciowy porz±dek *)
type porownanie = Mniejsze | Rowne | Wieksze;;
module type PORZADEK_LINIOWY = 
  sig
    type t
    val porownaj : t -> t -> porownanie
  end;;

(* Kolejka priorytetowa *)
module type PRI_QUEUE_FUNC = functor (Elem : PORZADEK_LINIOWY) -> 
  sig
    exception EmptyQueue
    type elem = Elem.t
    type queue
    val empty : queue 
    val is_empty : queue -> bool
    val put : elem -> queue -> queue
    val getmax : queue -> elem
    val removemax : queue -> queue
  end;;


module PriQueue : PRI_QUEUE_FUNC = functor (Elem : PORZADEK_LINIOWY) -> 
  struct
    type elem = Elem.t
    exception EmptyQueue
    type queue = elem list
    let empty = []
    let is_empty q = q = empty
    let rec put x q = 
      if q = [] then [x] 
      else if Elem.porownaj x (hd q) = Wieksze then x :: q
      else (hd q)::(put x (tl q))
    let getmax q =
      if q = [] then raise EmptyQueue
      else hd q
    let removemax (q:queue) = 
      if q = [] then raise EmptyQueue
      else tl q
  end;;

(* Kolejka priorytetowa liczb ca³kowitych *)
(* Uwaga: nie ograniczamy jej sygnatur± PORZADEK_LINIOWY ! *)
module IntOrder : (PORZADEK_LINIOWY with type t = int) = 
  struct
    type t = int
    let porownaj (x:int) (y:int) = 
      if x < y then Mniejsze else
      if x > y then Wieksze else Rowne
  end;;

module IntQueue = PriQueue (IntOrder);;

IntQueue.getmax (IntQueue.put 42 IntQueue.empty);;


(*
(* Poni¿sze definicje s± niepoprawne, gdy¿ zbyt du¿o zostaje przys³onione. *)
module SomeOrder : PORZADEK_LINIOWY = 
  struct
    type t = int
    let porownaj (x:int) (y:int) = 
      if x < y then Mniejsze else
      if x > y then Wieksze else Rowne
  end;;
module SomeQueue = PriQueue (SomeOrder);;
SomeQueue.getmax (SomeQueue.put 42 SomeQueue.empty);;
*)



module type S = sig val x : int end;;
module F (M : S) = struct let a = M.x let b = M.x end;;


(* Sortowanie *)
module type SORTING = functor (Elem : PORZADEK_LINIOWY) -> 
  sig
    type elem = Elem.t
    val sortuj : elem list -> elem list
  end;;

module HeapSort (PrQ : PRI_QUEUE_FUNC): SORTING = 
  functor (Elem : PORZADEK_LINIOWY) -> 
    struct
      type elem = Elem.t
      module PQ = PrQ (Elem)
      open PQ
      let wkladaj l =
        List.fold_left (fun h x -> put x h) empty l
      let rec wyjmuj q a = 
        if is_empty q then a 
        else wyjmuj (removemax q) (getmax q :: a)
      let sortuj l = wyjmuj (wkladaj l) []
    end;;

module Sortowanie = HeapSort (PriQueue);;
module S = Sortowanie (IntOrder);;
assert (S.sortuj [1;7;3;6;5;2;4;8] = [1; 2; 3; 4; 5; 6; 7; 8]);;

module Sort = HeapSort (PriQueue) (IntOrder);;
assert (Sort.sortuj [1;7;3;6;5;2;4;8] = [1; 2; 3; 4; 5; 6; 7; 8]);;

(* Wybór mediany *)
module type MEDIANA = 
  functor (Elem : PORZADEK_LINIOWY) -> 
    sig
      type elem = Elem.t
      exception PustaLista
      val mediana : elem list -> elem
    end;;

module Mediana (Sort : SORTING) : MEDIANA = 
    functor (Elem : PORZADEK_LINIOWY) -> 
      struct
        type elem = Elem.t
        exception PustaLista
        module S = Sort (Elem)
        let mediana l = 
          let s = S.sortuj l
          and n = List.length l
          in 
            if n = 0 then 
              raise PustaLista 
            else 
              List.nth s (n / 2)
      end;;

module M = Mediana (HeapSort (PriQueue)) (IntOrder);;
assert (M.mediana [4;3;5;6;2;1;7] = 4);;

(* Morfizmy sygnatur i redukty *)

module type S1 = 
  sig
    type t1
    type t2
    val a : t1
    val b : t2
    val f : t1 -> t2
  end;;

module type S2 =   
  sig
    type s1
    type s2
    val x : s1
    val y : s2
    val p : s1 -> s2
  end;;

module type S = 
  sig 
    type t
    val v : t
    val m : t -> t
  end;;

(* Wersja pierwsza -- nie widaæ co jest czym *)
module type M1 = functor (Sig : S) -> S1;;
module type M2 = functor (Sig : S) -> S2;;

module F1 : M1 = functor (Mod : S) ->
  struct
    open Mod
    type t1 = t
    type t2 = t
    let a = v
    let b = v
    let f = m
  end;;

module F2 : M2 = functor (Mod : S) ->
  struct
    open Mod
    type s1 = t
    type s2 = t
    let x = v
    let y = v
    let p = m
  end;;
    
module Suma = functor (M : S) -> 
  struct
    include F1(M)
    include F2(M)
  end;;

module Test = 
  Suma(struct
    type t = int
    let v = 42
    let m i = i
  end);;


(* Wersja druga -- widaæ co jest czym *)

module F1 = functor (Mod : S) ->
  struct
    open Mod
    type t1 = t
    type t2 = t
    let a = v
    let b = v
    let f = m
    let nowe1 = "42"
  end;;

module F2 = functor (Mod : S) ->
  struct
    open Mod
    type s1 = t
    type s2 = t
    let x = v
    let y = v
    let p = m
    let nowe2 = 42.0
  end;;
    
module Pushout = functor (M : S) -> 
  struct
    include F1(M)
    include F2(M)
  end;;

module Test = 
  Pushout(struct
    type t = int
    let v = 42
    let m i = i
  end);;


(* Przestrzenie metryczne *)
module type METRIC_SPACE = sig
  type t
  val d: t -> t -> float
end;;

module Manhattan_Plane : METRIC_SPACE = struct
  type t = float * float
  let d (x1, y1) (x2, y2) = abs_float (x1 -. x2) +. abs_float (y1 -. y2)
end;;

module Euclidian_Plane : METRIC_SPACE = struct
  type t = float * float
  let square x = x *. x
  let d (x1, y1) (x2, y2) = sqrt (square (x1 -. x2) +. square (y1 -. y2))
end;;

module Max_Plane : METRIC_SPACE = struct
  type t = float * float
  let d (x1, y1) (x2, y2) = max (abs_float (x1 -. x2)) (abs_float (y1 -. y2))
end;;

module Manhattan_Product (S1 : METRIC_SPACE) (S2 : METRIC_SPACE) : METRIC_SPACE = struct
  type t = S1.t * S2.t
  let d (x1, y1) (x2, y2) = S1.d x1 x2 +. S2.d y1 y2
end;;

module Euclidian_Product (S1 : METRIC_SPACE) (S2 : METRIC_SPACE) : METRIC_SPACE = struct
  type t = S1.t * S2.t
  let square x = x *. x
  let d (x1, y1) (x2, y2) = sqrt (square (S1.d x1 x2) +. square (S2.d y1 y2))
end;;

module Max_Product (S1 : METRIC_SPACE) (S2 : METRIC_SPACE) : METRIC_SPACE = struct
  type t = S1.t * S2.t
  let d (x1, y1) (x2, y2) = max (S1.d x1 x2) (S2.d y1 y2)
end;;

module Line : METRIC_SPACE = struct
  type t = float
  let d x y = abs_float (x -. y)
end;;

module Manhattan_Plane = Manhattan_Product (Line) (Line);;
module Euclidian_Plane = Euclidian_Product (Line) (Line);;
module Max_Plane       = Max_Product (Line) (Line);;

