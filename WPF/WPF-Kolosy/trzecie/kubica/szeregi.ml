(* Wersja z rekurencj± ogonow± *)
let szereg f n = 
  let rec sum a n = 
    if n = 0 then a else sum (a +. f n) (n - 1)
  in
    sum 0.0 n;;

(* Wersja bez rekurencji ogonowej *)
let rec szereg f n = 
  if n = 0 then 0.0 else f n +. szereg f (n-1);;


(* Przyk³adowy szereg *)

let szereg_pi_8 n = 
  szereg 
    (function i -> 
      1. /. ((4. *. float i -. 3.) *. (4. *. float i -. 1.)))
    n;;

let pi = 8. *. szereg_pi_8 1000;;
 

