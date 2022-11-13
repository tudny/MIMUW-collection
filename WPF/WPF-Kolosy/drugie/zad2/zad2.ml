

(* Rozwiązanie czas: O(n^2) | pamięć : O(n) *)
(* Dla każdego pasażera sprawdzam czy mógł się minąć z innym *)
(* Liczenie dla każdego zaczynam od -1, gdyż oczywiste jest, *)
(* że ten sposób policzy jakby pasażer A minął się z A, co nie jest prawdą *)
(* (nie jest prawdą tak długo, jak nie możemy minąć się sami) *)

let wirus li =
  let mijali (a, b) (c, d) =
    not (b <= c || d <= a) in
  List.fold_left (fun acc el ->
       List.fold_left (fun _acc el2 ->
           if mijali el el2 then _acc + 1
           else _acc
         ) (-1) li :: acc
    ) [] li |> List.rev
;;

assert ( wirus [(1,8);(2,5);(4,7);(3,4);(4,5);(6,7);(8,9)] = [5; 4; 4; 2; 3; 2; 0] );;

(*

type tree = Node of int * tree * tree | Null

let maks li =
  List.fold_left (fun _acc (a, b) -> max (max a b) _acc) 0 li
;;

let empty = Null
;;

let depth = 3;;

let rec empty lvls =
  if lvls = 0 then Null
  else Node (0, empty (lvls - 1), empty (lvls - 1))
;;

let get = function
  | Null -> 0
  | Node(v, l, r) -> v
;;

let dodaj (a, b) ile t =
  let rec add tr (akt_a, akt_b) =
    match tr with
    | Null -> Null
    | Node (va, l, p) ->
      if b < akt_a || akt_b < a then
        tr
      else if a <= akt_a && akt_b <= b then
        Node (va + ile, l, p)
      else
        let mid = (akt_a + akt_b) / 2 in
        let l = add l (akt_a, mid)
        and p = add p (mid + 1, akt_b)
        in Node(max (get l) (get p), l, p)
  in add t (1, 1 lsl (depth - 1))
;;

let czyt (a, b) t =
  let rec loop tr (akt_a, akt_b) =
    match tr with
    | Null -> 0
    | Node (va, l, p) ->
      if b < akt_a || akt_b < a then
        0
      else if a <= akt_a && akt_b <= b then
        va
      else
        let mid = (akt_a + akt_b) / 2 in
        let l = loop l (akt_a, mid)
        and p = loop p (mid + 1, akt_b)
        in max l p
  in loop t (1, 1 lsl (depth - 1))
;;

let t = empty depth;;

let t1 = dodaj (2, 4) 1 t;;

(*
let wirus li =
  let tr = List.fold_left (fun tr p -> dodaj p 1 tr) (empty depth) li in
  List.fold_left (fun _acc p -> czyt p tr :: _acc) [] li |> List.rev
;;


wirus [(1,8);(2,5);(4,7);(3,4);(4,5);(6,7);(8,9)];;


List.fold_left (fun tr p -> dodaj p 1 tr) (empty depth) [(1,8);(2,5);(4,7);(3,4);(4,5);(6,7);(8,9)];; *) *)
