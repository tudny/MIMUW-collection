(* Metoda Monte Carlo *)

let montecarlo dane wlasnosc liczba_prob = 
  let rec montecarlo_rec a n = 
    if n = 0 then 
      a
    else
      if wlasnosc (dane ()) then 
        montecarlo_rec (a+1) (n-1)
      else
        montecarlo_rec a (n-1)
  in
    float (montecarlo_rec 0 liczba_prob) /. float (liczba_prob);;
    
(* Generowanie liczb losowych *)

(* Wyszukiwanie du¿ych liczb pierwszych - na potrzeby sta³ych *)
let maxpierwsza n = 
  let pierwsza n = 
    let rec p x n = 
      if x * x > n then true
      else 
        if n mod x = 0 then false
        else
          p (x+1) n
    in 
      p 2 n
  in 
    let rec spr n = 
      if pierwsza n then 
        n 
      else 
        spr (n-1)
    in 
      spr n;;
      
(* Losowa liczba calkowita *)
let losowy_int = 
  let stan = ref 0
  and a = 937126433
  and b = 937187
  in function () -> 
    (
      stan := !stan * b + a;
      !stan
    );;


(* Losowa liczba rzeczywista nale¿±ca do przedzia³u [0, 1] *)
let losowa () = 
  let duzo = 1000
  in
    let l = float (losowy_int () mod (duzo + 1)) /. float (duzo)
    in 
      if l < 0.0 then (-. l) else l;;

(* Losowa liczba rzeczywista z podanego przedzia³u *)
let losowa_od_do od doo = losowa () *. (doo -. od) +. od;;


(* Zastosowanie metody Monte Carlo do przybli¿enia pi *)

(* Podnoszenie do kwadratu *)
let square x = x *. x;;

(* Ko³o jednostkowe *)
let kolo (x, y) = square x +. square y <= 1.0;;

(* Losowe dane nale¿±ce do kwadratu [-1, 1]x[-1, 1] *)
let kwadrat () = (losowa_od_do (-1.0) 1.0, losowa_od_do (-1.0) 1.0);;

let pi = 4.0 *. (montecarlo kwadrat kolo 1000);;

