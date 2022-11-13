(* Spamiêtywanie *)

#use "mapy.ml";;
open BstMap;;

(* Liczby Fibonacciego ze spamiêtywaniem, bez memoizatora. *)
let tab = ref empty;;
let rec fib n = 
  if dom !tab n then 
    apply !tab n
  else 
    let wynik = if n <= 1 then n else fib (n-1) + fib (n-2)
    in begin
      tab := update !tab n wynik;
      wynik
    end;;

let memoize tab f x = 
  if dom !tab x then 
    apply !tab x 
  else
    let wynik = f x 
    in begin
      tab := update !tab x wynik;
      wynik
    end;;

(* Liczby Fibonacciego ze spamiêtywaniem *)
let tab = ref empty;;
let rec fib n =
  memoize tab (function n -> 
    if n <= 1 then n else fib (n-1) + fib (n-2)) n;;


