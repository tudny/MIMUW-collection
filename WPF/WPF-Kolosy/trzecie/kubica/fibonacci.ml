(* Ró¿ne sposoby obliczania liczb Fibbonacciego. *)

(* Najprostsze, wyk³adnicze. *)
let rec fib n = 
  if n < 2 then n else fib (n - 1) + fib (n - 2);;

(* Na parach, bez def. lokalnych. *)
let rec fibpom n = 
  if n = 0 then 
    (0, 1)
  else
    let (a, b) = fibpom (n - 1)
    in (b, a + b);;

let fib n =
  let (a, b) = fibpom n
  in a;;

(* Na parach, z def. lokalnymi. *)
let fib n =
  let rec fibpom n = 
    if n = 0 then 
      (0, 1)
    else
      let (a, b) = fibpom (n - 1)
      in (b, a + b)
  in 
    let (a, b) = fibpom n
    in a;;

(* Pary w akumulatorach. *)
let fib n =
  let rec fibpom a b n = 
    if n = 0 then 
      a
    else
      fibpom b (a + b) (n - 1)
  in 
    fibpom 0 1 n;;

(* Logarytmicznie, przez potêgowanie macierzy *)
let parzyste n = n mod 2 = 0;;

let mnoz ((x11, x12), (x21, x22)) ((y11, y12), (y21, y22)) = 
  ((x11 * y11 + x12 * y21, x11 * y12 + x12 * y22), 
   (x21 * y11 + x22 * y21, x21 * y12 + x22 * y22));;

let square x = mnoz x x;;

let f = ((0, 1), (1, 1));;

let id = ((1, 0), (0, 1));;

let potega b n = 
  let rec pot b n a = 
    if n = 0 then a 
    else if parzyste n then pot (square b) (n / 2) a
    else pot b (n - 1) (mnoz a b)
  in
    pot b n id;;

let fib n = 
  let ((_, v), _) = potega f n 
  in v;;
      
    


