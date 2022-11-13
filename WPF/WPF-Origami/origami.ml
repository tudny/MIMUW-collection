(*******************************)
(**         ORIGAMI           **)
(**  Autor: Aleksander Tudruj **)
(** Review: Katarzyna Mielnik **)
(*******************************)


type point = float * float

type kartka = point -> int

(* Magiczny epsilon dokładności *)
let eps = 1e-9

(* Rzutowanie boola na inta co do wartości logicznej *)
let int_of_bool = function
  | true -> 1
  | false -> 0

(* Kwadrat liczby *)
let square x = x *. x

(* Funkcja porównująca dwie liczby typu float z dokładnościa do eps = 1e-9 *)
let compare a b =
  if abs_float (a -. b) <= eps then 0
  else if a < b then -1
  else 1

(* OPERATORY porównujące FLOATy *)

let (<.) a b =
  compare a b < 0

let (=.) a b =
  compare a b = 0

let (<=.) a b =
  a <. b || a =. b

(* Sprawdzenie czy liczba [x] jest w przedziale [] *)
let in_range x (a, b) =
  a <=. x && x <=. b

let prostokat (x1, y1 : point) (x2, y2 : point) (a, b : point) =
  int_of_bool
    (in_range a (x1, x2) &&
     in_range b (y1, y2))

let kolko (x, y : point) r (a, b : point) =
  let dx = (x -. a)
  and dy = (y -. b) in
  int_of_bool
    (square dx +. square dy <=. square r)

(* Wyznacznik macierzy 2x2 *)
let det (x, y) (a, b) =
  x *. b -. y *. a

(* Położenie punktu względem prostej *)
(* Na podstawie iloczynu wektorowego *)
let cmp (x1, y1) (x2, y2) (a, b) =
  det (x2 -. x1, y2 -. y1) (a -. x1, b -. y1)

(* Wyliczenie równania prostej, postaci ogólnej,
   przechodzącej przez dwa punkty *)
let get_line (x1, y1) (x2, y2) =
  let bB = x2 -. x1
  and aA = y1 -. y2
  and cC = x1 *. (y2 -. y1) -. y1 *. (x2 -. x1) in
  (aA, bB, cC)

(* Odbicie punktu względem prostej wyznaczonej przez dwa punkty *)
let reflect_point_over_line (a, b) p1 p2 =
  let (aA, bB, cC) = get_line p1 p2 in
    let sum = square aA +. square bB
    and diff = square aA -. square bB
    and dAB = 2. *. aA *. bB in
      let u_m = (-.diff) *. a -. dAB *. b -. 2. *. aA *. cC
      and v_m = diff *. b -. dAB *. a -. 2. *. bB *. cC in
        let u = u_m /. sum
        and v = v_m /. sum
  in (u, v)

let zloz p1 p2 (kar : kartka) q =
  match cmp p1 p2 q with
  | w when w  =. 0. -> kar q
  | w when w <=. 0. -> 0
  | _ ->
    let r = reflect_point_over_line q p1 p2 in
    kar q + kar r

let skladaj li kar =
  List.fold_left (fun acc (p1, p2) -> zloz p1 p2 acc) kar li
