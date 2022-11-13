(* Obliczanie pierwiastków metod± Newtona *)

(* Definicje ogólnie przydatnych rzeczy *)
let dodatnie x = x > 0.;;
let ujemne x = x < 0.;;
let square x = x *. x;;
let epsilon = 0.000001;;
let average x y = (x +. y) /. 2.0;;
let id x = x;;

(* Ró¿ñiczkowanie funkcji *)
let rozniczka f x dx = (f (x +. dx) -. f x) /. dx;;
let pochodna f x = rozniczka f x epsilon;;

(* Przypadek szczególny: pierwiastkowanie *)
let sqrt x = 
  let guess x = 1.0 
  and good g x = abs_float ((square g) -. x) <= epsilon *. g
  and improve g x = average g (x /. g) 
  in 
    let rec sqrt_iter g x = 
      if good g x then 
	g 
      else 
	sqrt_iter (improve g x) x
    in
      sqrt_iter (guess x) x;;


(* iteracyjny schemat poprawiania wyniku *)
let rec iteruj poczatek popraw czy_dobre wynik =
  if czy_dobre poczatek then 
    wynik poczatek
  else
    iteruj (popraw poczatek) popraw czy_dobre wynik;;


(* Szukanie zer funkcji przez bisekcjê *)
let bisekcja f l p = 
  let rec szukaj l p =
    let czy_dobre x = 
      abs_float (f (fst x)) < epsilon
    and popraw x = 
      let l = fst x
      and p = snd x
      in
        let srodek = average l p
	in
	  if dodatnie (f srodek) then 
	    (l, srodek) 
	  else
	    (srodek, p)
    in
      iteruj 
        (l, p)
        popraw
        czy_dobre
        fst
  in
    let wartosc_l = f l
    and wartosc_p = f p 
    in
      if ujemne wartosc_l && dodatnie wartosc_p then 
	szukaj l p
      else 
	szukaj p l;;

let sqrt a = 
  let f x = a -. square x
  in bisekcja f 0. (a +. 1.);;

(* Obliczanie \lfloor \sqrt n \rfloor przez bisekcjê wprost. *)
let square n = n * n;;
let sqrt n = 
  let rec bisekcja l p = 
    (* Niezmiennik: l <= \sqrt(n) < p+1 *)
    if l = p then l 
    else
      let sr = (l+p) / 2
      in 
        if square (sr+1) > n then 
	  bisekcja l sr
	else
	  bisekcja (sr+1) p
  in bisekcja 0 n;;

(* Szukanie zer metod± stycznych Newtona. *)
let newton f x = 
  let p = pochodna f
  in 
    let czy_dobre x = abs_float (f x) < epsilon
    and popraw x = x -. (f x) /. (p x)
    in
      iteruj x popraw czy_dobre id;;

(* Pierwiastkowanie metod± stycznych Newtona. *)
let sqrt a = 
  let f x = a -. square x
  in newton f 1.;;


(* Punkty sta³e funkcji. *)
let punkt_staly f x = 
  let blisko x = abs_float (x -. f x) < epsilon
  in
    iteruj x f blisko f;;

(* Wyt³umianie przez u¶rednianie. *)
let usrednienie f = 
  function x -> average x (f x);;

let sqrt x = 
  punkt_staly (usrednienie (function y -> x /. y)) 1.0;;
