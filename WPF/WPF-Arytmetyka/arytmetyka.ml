(* Autor Aleksander Tudruj *)
(* Code review Kacper Jabłoński *)

(* typ 'wartosc' opisana w treści zadania *)
(* zawiera dwa konstruktory *)
(* Przedzial (a, b) - przechowuje informacje o zbiorze [a, b] *)
(* Dopelnienie (x, y) - przechowuje informacje o zbiorze [-inf, x] u [y, inf] *)

type wartosc = Przedzial of float * float |
               Dopelnienie of float * float |
               Pusty;;

let is_nan x =
  compare x nan = 0

(* ****************** *)
(*    KONSTRUKTORY    *)
(* ****************** *)

let wartosc_dokladnosc x p =
  let dokladnosc = abs_float (x *. p /. 100.) in
  Przedzial (x -. dokladnosc, x +. dokladnosc)

let wartosc_od_do x y =
  Przedzial (x, y)

let wartosc_dokladna x =
  wartosc_dokladnosc x 0.

(* *************** *)
(*    SELEKTORY    *)
(* *************** *)

let in_wartosc w x =
  match w with
  | Przedzial (a, b)   -> (a <= x && x <= b)
  | Dopelnienie (a, b) -> (x <= a || b <= x)
  | Pusty              -> false

let min_wartosc w =
  match w with
  | Przedzial (a, _)   -> a
  | Dopelnienie (_, _) -> neg_infinity
  | Pusty              -> nan

let max_wartosc w =
  match w with
  | Przedzial (_, b)   -> b
  | Dopelnienie (_, _) -> infinity
  | Pusty              -> nan

let sr_wartosc w =
  ((min_wartosc w) +. (max_wartosc w)) /. 2.

(* ****************** *)
(*    MODYFIKATORY    *)
(* ****************** *)

(* DEBUG *)
let rec print_list_of_floats li =
  match li with
  | [] -> print_endline ""
  | (h :: t) -> print_float h; print_string " "; print_list_of_floats t;;

(* DEBUG *)
let print_wartosc w =
  match w with
  | Przedzial (a, b) -> print_string "Przedzial: "; print_float a; print_string " "; print_float b; print_endline "";
  | Dopelnienie (a, b) -> print_string "Dopelnienie: "; print_float a; print_string " "; print_float b; print_endline ""
  | Pusty -> print_string "Pusty"; print_endline ""

(* DEBUG *)
let print_pair_of_floats (a, b) =
  print_string "("; print_float a; print_string ", "; print_float b; print_string ")\n"

(* Funkcja upraszczająca utworzone przedziały, bądź doprowadzająca je do poprawności. Jeżeli nie znajdzie nic do uproszczenia zwraca przekazany przedział *)
let normalize w =
  match w with
  | Dopelnienie (a, b) when a >= b           -> Przedzial (neg_infinity, infinity)
  | Dopelnienie (a, b) when a = neg_infinity -> Przedzial (b, infinity)
  | Dopelnienie (a, b) when b = infinity     -> Przedzial (neg_infinity, a)
  | Dopelnienie (a, b) when a = neg_infinity && b = infinity -> Pusty
  | Przedzial (a, b)   when a = nan || b = nan               -> Pusty
  | _ -> w

(* Komparator dwóch floatów (x y) za pomocą porównania f (float -> float -> bool) *)
(* Wartość nan zwracana jest tylko gdy x = nan = y *)
let num_comp f x y =
  if is_nan x then y else if is_nan y then x else
  if f x y then x
  else y;;

(* Minimum floatów *)
let min_num =
  num_comp (<);;

(* Maksimum floatów *)
let max_num =
  num_comp (>);;

(* Element ekstremalny na liście (np min lub max) *)
let rec extremum_in_list comparator li =
  match li with
  | []  -> nan
  | [x] -> x
  | h :: t -> comparator h (extremum_in_list comparator t);;

(* Minium z listy *)
let min_in_list =
  extremum_in_list (min_num)

(* Maksimum z listy *)
let max_in_list =
  extremum_in_list (max_num)

let rec plus x y =
  match x, y with
  | Przedzial (a, b), Przedzial (c, d) ->
    let kandydaci = (a +. c) :: (a +. d) :: (b +. c) :: (b +. d) :: [] in
    normalize (Przedzial (min_in_list kandydaci, max_in_list kandydaci))
  | Dopelnienie (_, _), Dopelnienie (_, _) -> Przedzial (neg_infinity, infinity)
  | Dopelnienie (_, _), Przedzial (_, _) -> plus y x
  | Przedzial (a, b), Dopelnienie (c, d) ->
      normalize (Dopelnienie (max_num (a +. c) (b +. c), min_num (a +. d) (b +. d)))
  | _, Pusty -> Pusty
  | Pusty, _ -> Pusty

(* Dla danego zbioru "w" zwraca -w = {-x | x \in w} *)
let przeciwna w =
  match w with
    | Przedzial (a, b) -> normalize (Przedzial (~-.b, ~-.a))
    | Dopelnienie (a, b) -> normalize (Dopelnienie (~-.b, ~-.a))
    | Pusty -> Pusty

let minus x y = plus x (przeciwna y);;

(* Zwraca znak liczby *)
let znak a =
  match classify_float a with
  | FP_normal -> a /. abs_float(a)
  | FP_subnormal -> a /. abs_float(a)
  | FP_zero -> 0.
  | FP_infinite -> if a = infinity then 1. else -.1.
  | FP_nan -> nan

(* sprawdza czy a i b są po różnych stronach osi Oy *)
let czy_przeciwne_znaki a b =
  ((znak a) *. (znak b)) < 0.;;

let czy_zawiera_zero_nie_na_krancu w =
  match w with
  | Przedzial (a, b) -> czy_przeciwne_znaki a b
  | Dopelnienie (a, b) -> not (czy_przeciwne_znaki a b)
  | Pusty -> false

let rec razy x y =
  let zero = wartosc_dokladna 0. in
  match x, y with
  | Pusty, _ -> Pusty
  | _, Pusty -> Pusty
  | x, y when x = zero || y = zero -> zero
  | Przedzial (a, b), Przedzial (c, d) ->
      let kandydaci = (a *. c) :: (a *. d) :: (b *. c) :: (b *. d) :: [] in
      normalize (Przedzial (min_in_list kandydaci, max_in_list kandydaci))
  | Przedzial (a, b), Dopelnienie (c, d) ->
      let w1 = razy x (Przedzial (neg_infinity, c))
      and w2 = razy x (Przedzial (d, infinity)) in
      let (mini1, maks1) = (min_wartosc w1, max_wartosc w1)
      and (mini2, maks2) = (min_wartosc w2, max_wartosc w2) in
      let (mini, maks) = if maks1 = infinity then (maks2, mini1) else (maks1, mini2) in
      normalize (Dopelnienie (mini, maks))
  | Dopelnienie (_, _), Przedzial (_, _) -> razy y x
  | Dopelnienie (a, b), Dopelnienie (c, d) ->
      let w1 = razy (normalize (Przedzial (neg_infinity, a))) y
      and w2 = razy (normalize (Przedzial (b, infinity))) y in
      match w1, w2 with
      | Przedzial (_, _), _ -> w1
      | _, Przedzial (_, _) -> w2
      | Dopelnienie (a, b), Dopelnienie (c, d) -> normalize (Dopelnienie (max_num a c, min_num b d))
      | _ -> Pusty

(* Dla danego zbioru "w" zwraca w^{-1} = {x^{-1} | x \in w\{0}} *)
let odwrotnosc w =
  match w with
  | Pusty -> Pusty
  | Przedzial (a, b) when a = 0. && b = 0. -> Pusty
  | Przedzial (a, b) when a = 0. -> normalize (Przedzial (1. /. b, infinity))
  | Przedzial (a, b) when b = 0. -> normalize (Przedzial (neg_infinity, 1. /. a))
  | Przedzial (a, b) when czy_zawiera_zero_nie_na_krancu w -> normalize (Dopelnienie (1. /. a, 1. /. b))
  | Przedzial (a, b) -> normalize (Przedzial (1. /. b, 1. /. a))
  | Dopelnienie (a, b) when a = 0. && b = 0. -> w
  | Dopelnienie (a, b) when a = 0. -> normalize (Przedzial (neg_infinity, 1. /. b))
  | Dopelnienie (a, b) when b = 0. -> normalize (Przedzial (1. /. a, infinity))
  | Dopelnienie (a, b) when czy_zawiera_zero_nie_na_krancu w -> normalize (Dopelnienie (1. /. b, 1. /. a))
  | Dopelnienie (a, b) -> normalize (Przedzial (1. /. a, 1. /. b))

let podzielic x y =
  razy x (odwrotnosc y);;
