(* Stopieñ parzysto¶ci danej liczby. *)

(* Proste rozwi±zanie nieogonowe. *)
(* T(n) = M(n) = O(log n)         *)
let rec stopien_parzystosci n = 
  if n = 0 then -1 
  else if n mod 2 = 0 then 1 + stopien_parzystosci (n/2) 
  else 0;;

assert (List.map stopien_parzystosci [0;1;2;3;4;5;6;12;24;32]=[-1; 0; 1; 0; 2; 0; 1; 2; 3; 5]);;	


(* Proste rozwi±zanie ogonowe.  *)
(* T(n) = O(log n), M(n) = O(1) *)
let stopien_parzystosci n = 
  let rec pom a n = 
    if n mod 2 = 0 then pom (a+1) (n/2)
    else a
  in
    if n = 0 then -1 else pom 0 n;;

assert (List.map stopien_parzystosci [0;1;2;3;4;5;6;12;24;32]=[-1; 0; 1; 0; 2; 0; 1; 2; 3; 5]);;	


(* Rozwi±zanie szybsze. *)
(* T(n) = O((log log n)^2), M(n) = O(1) *)

let stopien_parzystosci n = 
  (* p = 2^k | n  lub n nieparzyste,  pom a k n = a + stopien_parzystosci n *)
  let rec pom a k p n = 
    if n mod 2 <> 0 then a 
    else 
      let s = p*p 
      in
        if n mod s = 0 then pom a (2*k) s n 
	else  (assert (n mod p = 0); pom (a+k) 1 2 (n / p))
  in 
    if n = 0 then -1 else pom 0 1 2 n;;

assert (List.map stopien_parzystosci [0;1;2;3;4;5;6;12;24;32]=[-1; 0; 1; 0; 2; 0; 1; 2; 3; 5]);;	

(* Rozwi±zanie optymalne. *)
(* T(n) = M(n) = O(log log n) *)
let stopien_parzystosci n = 
  (* l = [2^k; 2^(k-1); ...; 2; 1], pom a k l n = a + stopien_parzystosci n. *)
  let rec pom a k ((h::t) as l) n = 
    if n mod 2 <> 0 then a 
    else 
      let s = h * h 
      in
        if n mod s = 0 then pom a (2*k) (s::l) n 
	else if n mod h = 0 then pom (a+k) (k/2) t (n/h)
	else pom a (k/2) t n
  in 
    if n = 0 then -1 else pom 0 1 [2; 1] n;;

assert (List.map stopien_parzystosci [0; 1; 2; 3; 4; 5; 6; 12; 24; 32; 896]=[-1; 0; 1; 0; 2; 0; 1; 2; 3; 5; 7]);;	
