(* Algorytm Archimedesa przybliżania liczby π.                          *)
(* Rozważamy n-kąt foremny wpisany i opisany na okręgu jednostkowym.    *)
(* Długość boku n-kąta wpisanego oznaczmy przez a, a opisanego przez b. *)

open List;;

(* Wykonuje k kroków algorytmu Archimedesa, zaczynając od 6-kąta, a kończąc na 6*2^k kącie. *)
(* Wynikiem jest lista kolejno uzyskanych krotek: (n, a, b, dolne prz. pi, górne prz. pi)   *)
let archimedes k = 
  let rec iter k n a b acc = 
    if k = 0 then rev acc
    else 
      let a1 = sqrt(2. -. 2.*.a /. b)
      and b1 = 2. *. sqrt((b -. a) /. (b +. a))
      in iter (k-1) (2*n) a1 b1 ((2*n, a1, b1, float(n) *. a1, float(n) *. b1)::acc)
  in
    iter k 6 1. (2. /. sqrt 3.) [];;
