(* Prostsze algorytmy sprawdzaj±ce, czy liczba jest pierwsza. *)

(* Przydatne drobiazgi *)
let square x = x * x;;
let parzyste x = x mod 2 = 0;;

(* Prosty algorytm O(\sqrt(n)). *)
let min_dzielnik n =
  let rec dziel k =
    if square k > n then n
    else if (n mod k) = 0 then k
    else dziel (k + 1)
  in dziel 2;;



(* Test Millera-Rabina na bycie liczb± pierwsz±. *)

let min_dzielnik n =
  let rec dziel k =
    if square k > n then n
    else if (n mod k) = 0 then k
    else dziel (k + 1)
  in dziel 2;;


exception Pierwiastek;;

(* Potêgowanie modulo n z wykrywaniem nietrywialnych pierwiastków *)
let rec expmod b k n = 
  let test x = 
    if not (x = 1) && not (x = n - 1) && (square x) mod n = 1 then 
      raise Pierwiastek
    else
      x
  in
    if k = 0 then 1 else 
    if parzyste k then (square (test (expmod b (k / 2) n))) mod n
    else ((test (expmod b (k-1) n)) * b) mod n;;

(* Jednokrotny test losowy na bycie liczb± pierwsz±. *)
let randtest n = 
  if parzyste n then 
    n = 2 
  else if n = 3 then true
  else
    try 
      expmod (Random.int (n-3) + 2) (n-1) n = 1
    with Pierwiastek -> false;;
