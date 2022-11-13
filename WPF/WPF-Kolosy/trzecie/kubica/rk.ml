(* Algorytm Rabina-Karpa *)

open Array;;
open Char;;

let is_prime p = 
  let rec pom d = 
    if d * d > p then true
    else if p mod d = 0 then false 
    else pom (d+1)
  in pom 2;;

let rec primes n = 
  if is_prime n then begin
    print_int n;
    print_newline()
  end;
  if n > 2 then primes (n-1);;

let p = 257
and q = 4177969;;

(* Potęgowanie b^k mod q *)
let exp_mod b k = 
  let rec iter k a = 
    if k = 0 then a else iter (k - 1) ((a * b) mod q)
  in
    iter k 1;;

let fingerprint str k = 
  let rec sumuj acc pow pos = 
    if pos < 0 then acc 
    else sumuj ((acc + pow * code str.(pos)) mod q) ((pow * p) mod q) (pos -1)
  in sumuj 0 1 (k - 1);;

let shift fi pot x y = 
  let res = (((fi * p) mod q + y - x * pot) mod q)
  in
    if res < 0 then (res + p * q) mod q
    else res;;

let occurrence w t i = 
  let rec check j = 
    if j < 0 then true
    else if w.(j) = t.(i+j) then 
      check (j-1)
    else false
  in check (Array.length w - 1);;

let pattern_match w t = 
  let n = Array.length t
  and m = Array.length w
  in
    let f = fingerprint w m
    and pot = exp_mod p m 
    in
      let rec iter acc i fi = 
	let a = if (fi = f) && (occurrence w t i) then i::acc else acc 
	in 
	  if i = n-m then List.rev a
	  else iter a (i+1) (shift fi pot (code t.(i)) (code t.(i + m)))
      in iter [] 0 (fingerprint t m);;

assert (pattern_match [|'a'; 'l'; 'a'|] [|'b'; 'a'; 'l'; 'a'; 'l'; 'a'; 'n'; 'g'; 'a'; 'l'; 'a'|] = [1;3;8]);;
