open Lazy;;

(* Implementacja struktur danych za pomoc± procedur *)

(* Warto¶ci logiczne *)
let prawda x y = x;;
let falsz  x y = y;;
let jesli w k a = force (w k a);;
(* jesli prawda (lazy 1) (lazy 2);;  *)

let i x y  = x y falsz;;
let lub x y = x prawda y;;
let nie x a b = x b a;; 
(* jesli (lub falsz prawda) (lazy 1) (lazy 2);; *)

(* Produkty kartezjañskie *)
let pr x y = function f -> f x y;;
let p1 p = p prawda;;
let p2 p = p falsz;;

(* Parametr x jest potrzebny ze wzglêdu na s³abe zmienne typowe. *)
let p x = pr 1 "ala" x;;
p1 p;;
p2 p;;

(* Listy nieskoñczone *)
let stala x n =  x;;
let sklej h t n = if n = 0 then h else t (n-1);; 
let glowa l = l 0;;
let ogon l n = l (n + 1);;

(* Liczby naturalne Churcha *)
let id x = x;;
let compose f g = function x -> f (g x);;
let zero f = id;;
let jeden f = f;;

let ile n = n (function x -> x + 1) 0;;
let rec iterate n f =
  if n = 0 then id else compose (iterate (n-1) f) f;;

assert (ile zero = 0);;
assert (ile jeden = 1);;
assert (ile (iterate 1000) = 1000);;

let inc n f = compose (n f) f;;
let plus m n f = compose (m f) (n f);;
let razy = compose;;
let potega m n = n m;;

let dwa f = inc jeden f;;
let trzy f = plus dwa jeden f;;
let szesc f = razy dwa trzy f;;
let siedem f = inc szesc f;;
let osiem f = potega dwa trzy f;;

assert (ile (razy dwa dwa) = 4);;
assert (ile (razy szesc siedem) = 42);;
assert (ile (plus dwa (razy (plus dwa trzy) osiem)) = 42);;

let test_zera x = x (function _ -> falsz) prawda;;

(* jesli (test_zera zero) (lazy 1) (lazy 2);; *)

let dec n f x = 
  let (y, _) = n (function (a, b) -> (b, f b)) (x, x) 
  in y;;

let dec n f x = 
  p1 (n (function p -> pr (p2 p) (f (p2 p))) (pr x x));;

(* ile (dec dwa);; *)

let minus m n = (n dec) m;;

(* ile (minus dwa jeden);; *)

let wieksze m n = nie (test_zera (minus m n));;

(* jesli (wieksze dwa jeden) (lazy 1) (lazy 2);; *)
