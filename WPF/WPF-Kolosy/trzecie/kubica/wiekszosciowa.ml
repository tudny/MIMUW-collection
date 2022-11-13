(* Zadanie o liście większościowej. *)

open List;;

let majority l = 
  let rec scan c k l = 
    match l with
      [] -> c |
      h::t -> 
	if k = 0 then scan h 1 t 
	else if c = h then scan c (k+1) t
	else scan c (k-1) t
  in scan (hd l) 0 l;;
