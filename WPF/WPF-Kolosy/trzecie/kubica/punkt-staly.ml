(* Punkt sta�y rekurencji *)

(* Zwyk�a rekurencyjna silnia *)
let rec fac n = 
  if n <= 1 then 
    1 
  else 
    n * fac (n - 1);;

(* Definicja silni bez rekurencji. *)  
let factorial fac n = 
  if n <= 1 then 
    1 
  else 
    fac (n - 1);;

(* Operator punktu sta�ego -- pierwsze podej�cie *)
let rec y f = f (y f);;

let f = y factorial ;;
(* To sie zapetli *)

(* Podej�cie drugie, z uleniwieniem fac *)
let factorial fac n = 
  if n <= 1 then 
    1 
  else 
    n * (force fac (n - 1));;

let rec y f = f (lazy (y f));;

let fac = y factorial ;;

fac 3;;
