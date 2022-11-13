(******************************)
(* Autor: Aleksander Tudruj   *)
(* Review: Patryk Jędrzejczak *)
(******************************)


(* Typ wariantowy *)
(* Null -> Puste drzewo o prawej wysokości 0 *)
(* Node(l, v, r, h) -> Niepuste drzewo o prawej wysokości [h],
   wartości w wierzchołku [v] oraz podrzewach [l] lewym i [r] prawym *)

type 'a queue = Node of ('a queue * 'a * 'a queue * int) | Null;;

(* Zwraca pustą koleję *)
let empty = Null;;

(* Zwraca kolejkę (drzewo) z jednym elementem [x]
   z zachowaniem warunków prawej wysokości *)
let singleton x = Node (Null, x, Null, 1);;

(* Wyciąga wartość [prawa wysokość] z dowolnego drzewa *)
let get_height t =
  match t with
  | Null -> 0
  | Node(_, _, _, h) -> h

(* Łączy dwie kolejki *)
let rec join q1 q2 =
  match q1, q2 with
  | Null, Null -> Null
  | Null, _ -> q2
  | _, Null -> q1
  | Node(_, v1, _, _), Node(_, v2, _, _) when v1 > v2 -> join q2 q1
  | Node(l1, v1, r1, h1), Node(l2, v2, r2, h2) ->
    let q3 = join r1 q2 in
    let h_l1 = get_height l1
    and h_q3 = get_height q3 in
    if h_q3 > h_l1 then
      Node(q3, v1, l1, h_l1 + 1)
    else
      Node(l1, v1, q3, h_q3 + 1)

(* Dodale element [e] do kolejki [q] za pomocą złączenia drzewa
   z jednym elementem (e) z q*)
let add e q =
  join (singleton e) q

(* Wyjątek zwracany, gdy kolejka jest pusta *)
exception Empty

(* Zwraca najmniejszy element wraz z kolejką bez tego elementu *)
let delete_min = function
  | Null -> raise Empty
  | Node(l, v, r, _) -> (v, join l r)

(* Sprawdza czy kolejka jest pusta *)
let is_empty = function
  | Null -> true
  | Node(_) -> false
