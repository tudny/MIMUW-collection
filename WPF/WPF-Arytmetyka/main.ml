
open Arytmetyka;;

let is_nan a = compare a nan = 0;;

let zero = wartosc_dokladna 0.;;
let jeden = wartosc_dokladna 1.;;
let prosty = wartosc_od_do 0. 1.;;

let a = podzielic prosty zero;; (* <0, 1> / <0, 0> = Pusty *)
let b = podzielic prosty jeden;; (* <0, 1> / <1, 1> = <0, 1> *)
let c = podzielic zero prosty;; (* <0, 0> / <0, 1> = <0, 0> *)
let d = podzielic jeden prosty;; (* <1, 1> / <0, 1> = <1, inf> *)
let e = podzielic zero jeden;; (* <0, 0> / <1, 1> = <0, 0> *)
let f = podzielic jeden zero;; (* <1, 1> / <0, 0> = Pusty *)

assert ( (is_nan (max_wartosc a), is_nan(min_wartosc a), is_nan(sr_wartosc a)) = (true, true, true) );
assert ( (max_wartosc b, min_wartosc b, sr_wartosc b) = (1., 0., 0.5) );
assert ( (max_wartosc c, min_wartosc c, sr_wartosc c) = (0., 0., 0.) );
assert ( (max_wartosc d, min_wartosc d, sr_wartosc d) = (infinity, 1., infinity) );
assert ( (max_wartosc e, min_wartosc e, sr_wartosc e) = (0., 0., 0.) );
assert ( (is_nan (max_wartosc f), is_nan(min_wartosc f), is_nan(sr_wartosc f)) = (true, true, true) );
