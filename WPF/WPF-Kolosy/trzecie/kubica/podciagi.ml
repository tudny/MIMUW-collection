(* Znajdowanie najwiêkszego wspólnego podci±gu dwóch ci±gów *)
open Array;;
open List;;
#use "listy.ml";;

(************************************************************)
(* Rozwi±zanie imperatywne *)


(* Najd³u¿szy wspólny podci±g dwóch list. *)
let nwp ciag_x ciag_y = 
  if (ciag_x = []) || (ciag_y = []) then [] else
  let n = length ciag_x 
  and m = length ciag_y 
  and x = of_list ((hd ciag_x)::ciag_x)
  and y = of_list ((hd ciag_y)::ciag_y)
  in
    (* a.(i).(j) = lengt (nwp [x_1;..;x_i] [y_1;..;y_j]) *)
    let a = make_matrix (n+1) (m+1) 0 
    in 
      (* Rekonstrukcja nwp na podstawie tablicy d³ugo¶ci. *)
      let rec rekonstrukcja acc i j = 
	if i = 0 or j = 0 then acc
	else if x.(i) = y.(j) then 
	  rekonstrukcja (x.(i)::acc) (i-1) (j-1)
	else if a.(i).(j) = a.(i-1).(j) then 
	  rekonstrukcja acc (i-1) j
	else 
	  rekonstrukcja acc i (j-1)
      in
        begin
          (* Wype³nienie tablicy a *)
          for i = 1 to n do 
	    for j = 1 to m do 
	      if x.(i) = y.(j) then 
	        (* Odpowiadaj±ce sobie elementy ci±gów s± równe. *)
	        a.(i).(j) <- a.(i-1).(j-1) + 1
	      else
	        (* Wybór wariantu realizuj±cego d³u¿szy podci±g. *)
	        a.(i).(j) <- max a.(i-1).(j) a.(i).(j-1)
	    done
          done;
          (* Rekonstrukcja wyniku. *)
          rekonstrukcja [] n m 
          
        end;;

assert
  (let l = nwp ['A'; 'B'; 'C'; 'B'; 'D'; 'A'; 'B'] 
              ['B'; 'D'; 'C'; 'A'; 'B'; 'A']     
  in 
    (l = ['B'; 'D'; 'A'; 'B']) || 
    (l = ['B'; 'C'; 'A'; 'B']) ||
    (l = ['B'; 'C'; 'B'; 'A']));;


(************************************************************)
(* Rozwiazanie z mapami *)
#use "mapy.ml";;
open BstMap;;

(* Wzór na zale¿no¶æ dlugo¶ci w tablicy  *)
let dlugosc u v w x y = 
  if x = y then 
    v + 1
  else 
    max u w;;
    
(* 
  Zbudowanie tablicy wspólnych podci±gów 
  A (i, j) = d³ugo¶æ najd³u¿szego wspólnego podci±gu 
  (x_1, ..., x_i) i (y_1, ..., y_j). 
  Tablica jest uzupe³niona stra¿nikami równymi 0.
*)
let zbuduj_tablice ciag_x ciag_y = 

  (* 
     Dodaje kolejny wiersz tablicy, 
     y= rozpatrywany element ci±gu y, j = nr wiersza.
  *)
  let dodaj_wiersz tablica y j =
    let (t, _) = 
      fold_left
	(fun (tab, i) x ->
          let u = apply tab (j, i-1)
          and v = apply tab (j-1, i-1)
          and w = apply tab (j-1, i)
          in
            (update tab (j, i) (dlugosc u v w x y), i+1))
	(update tablica (j, 0) 0, 1)
	ciag_x
    in t
  in 
    (* Pierwszy wiersz tablicy -- same zera *)
    let pierwszy = 
      let (a, _) = 
	fold_left
	  (fun (a, i) _ -> (update a (0, i) 0, i+1))
	  (update empty (0, 0) 0, 1)
	  ciag_x
      in a
        
    (* Budowanie tablicy, yl= pozostale elementy ci±gu y, j = nr wiersza. *)
    and buduj a yl = 
      let (w, _) = 
	fold_left (fun (a, j) y -> (dodaj_wiersz a y j, j+1)) (a, 1) yl
      in w
    in
      buduj pierwszy ciag_y;;
    
    
(* Rekonstrukcja ci±gu na podstawie tablicy *)
let odtworz_podciag tablica ciag_x ciag_y = 
  let rec odtworz akumulator tab ciag_x ciag_y i j = 
    match ciag_x with 
      [] -> akumulator |
      (x::tx) -> 
        match ciag_y with 
          [] -> akumulator |
          (y::ty) -> 
            if x = y then 
              odtworz (x::akumulator) tab tx ty (i-1) (j-1) 
            else 
              let u = apply tab (j, i-1) 
              and w = apply tab (j-1, i)
              in 
                if u > w then 
                  odtworz akumulator tab tx ciag_y (i-1) j
                else 
                  odtworz akumulator tab ciag_x ty i (j-1)
  in 
    odtworz [] tablica 
            (rev ciag_x) (rev ciag_y)
            (length ciag_x) (length ciag_y);;
    
(* Najwiêkszy wspólny podci±g *)
let nwp ciag_x ciag_y = 
  let 
    a = zbuduj_tablice ciag_x ciag_y 
  in 
    odtworz_podciag a ciag_x ciag_y;;

assert
  (let l = nwp ['A'; 'B'; 'C'; 'B'; 'D'; 'A'; 'B'] 
              ['B'; 'D'; 'C'; 'A'; 'B'; 'A']     
  in 
    (l = ['B'; 'D'; 'A'; 'B']) || 
    (l = ['B'; 'C'; 'A'; 'B']) ||
    (l = ['B'; 'C'; 'B'; 'A']));;


(******************************************************************)

(* Implementacja za pomoc± struktury wska¼nikowej *)
open List;;
#use "listy.ml";;

(* Struktura tablicy *)
type array = { value : int; up : array; upleft : array; left : array };;

(* Same zera *)
let rec zero = { value = 0; up = zero; upleft = zero; left = zero };;

(* Wzór na zale¿no¶æ dlugo¶ci w tablicy  *)
let dlugosc u v w x y = 
  if x = y then 
    v + 1
  else 
    max u w;;
    
(* 
  Zbudowanie tablicy wspólnych podci±gów 
  A[i, j] = d³ugo¶æ najd³u¿szego wspólnego podci±gu 
  (x_1, ..., x_i) i (y_1, ..., y_j). 
  Tablica jest reprezentowana w postaci ODWROCONEJ listy ODWROCONYCH list.
  Tablica jest uzupe³niona stra¿nikami równymi 0.
*)
let zbuduj_tablice ciag_x ciag_y = 

  (* Dodaje kolejny wiersz tablicy, 
     y= rozpatrywany element ci±gu y.
  *)
  let dodaj_wiersz tab y  =
    let rec dodawaj tab_up lista = 
      match lista with
        []     -> zero |
        (x::t) -> 
          let tab_upleft = tab_up.left
          in 
            let tab_left = dodawaj tab_upleft t
            in 
	      let v = dlugosc tab_left.value tab_upleft.value tab_up.value x y
	      in
                { value  = v;
		  left   = tab_left;
		  upleft = tab_upleft;
		  up     = tab_up 
		}
    in 
      dodawaj tab (rev ciag_x)
  in 
    fold_left dodaj_wiersz zero ciag_y;;
    
    
(* Rekonstrukcja ci±gu na podstawie tablicy *)
let odtworz_podciag tablica ciag_x ciag_y = 
  let rec odtworz akumulator tab ciag_x ciag_y = 
    match ciag_x with 
      [] -> akumulator |
      (x::tx) -> 
        match ciag_y with 
          [] -> akumulator |
          (y::ty) -> 
            if x = y then 
              odtworz (x::akumulator) tab.upleft tx ty 
            else if tab.left.value > tab.up.value then 
              odtworz akumulator tab.left tx ciag_y
            else 
              odtworz akumulator tab.up ciag_x ty
  in 
    odtworz [] tablica (rev ciag_x) (rev ciag_y);;
    
(* Najwiêkszy wspólny podci±g *)
let nwp ciag_x ciag_y = 
  let 
    a = zbuduj_tablice ciag_x ciag_y 
  in 
    odtworz_podciag a ciag_x ciag_y;;


assert
  (let l = nwp ['A'; 'B'; 'C'; 'B'; 'D'; 'A'; 'B'] 
              ['B'; 'D'; 'C'; 'A'; 'B'; 'A']     
  in 
    (l = ['B'; 'D'; 'A'; 'B']) || 
    (l = ['B'; 'C'; 'A'; 'B']) ||
    (l = ['B'; 'C'; 'B'; 'A']));;

(***************************************************************)

(* Rozwi±zanie ze spamiêtywaniem.  *)
open List;;
#use "spamietywanie.ml";;

let nwp ciag_x ciag_y = 
  let x = Array.of_list ciag_x
  and y = Array.of_list ciag_y
  and tab = ref empty
  in
    let rec pom (i, j) = 
      memoize tab (fun (i, j) -> 
	if i = 0 or j = 0 then ([], 0)
	else if x.(i-1) = y.(j-1) then 
	  let (c, l) = pom (i-1, j-1)
	  in (x.(i-1)::c, l+1)
	else
	  let (c1, l1) = pom (i-1, j)
	  and (c2, l2) = pom (i, j-1)
	  in
	    if l1 > l2 then (c1, l1) else (c2, l2)
      ) (i, j)
    in
      rev (fst (pom (length ciag_x, length ciag_y)));;

assert
  (let l = nwp ['A'; 'B'; 'C'; 'B'; 'D'; 'A'; 'B'] 
              ['B'; 'D'; 'C'; 'A'; 'B'; 'A']     
  in 
    (l = ['B'; 'D'; 'A'; 'B']) || 
    (l = ['B'; 'C'; 'A'; 'B']) ||
    (l = ['B'; 'C'; 'B'; 'A']));;

(*************************************************************************)



