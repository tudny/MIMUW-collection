(* Algorytm mini-max z alfa-beta obci�ciem. *)

(* Gracze: max = 1, min = -1. *)
let i = 1;;

(* Przeciwnik gracza who *)
let other who = - who;;

(* Wyniki gier. *)
let won = 1;;
let lost = -1;;

(* Gra w kamienie - cz�� zale�na od gry *)
#use "gra-w-kamienie.ml";;

(* Algorytm alfa-beta obci�cia.                                   *)

(* Ograniczenie dolne wyniku gracza who *)
let no_result who = -2 * who;;

(* Obci�cie jest realizowane jako wyj�tek *)
exception Obciecie of int * move;;

(* Ocena sytuacji s, gdy pierwszy ruch wykonuje who, a przeciwnik *)
(* jak dot�d uzyska� wag� alpha.                                  *)
(* Wynikiem jest para (waga, ruch realizuj�cy wag�).              *)
let rec alfabeta s who alpha = 
  if finished s who then 
    (* Gra sko�czona *)
    (result s who, pass)
  else
    (* Gra trwa *)

    (* Funkcja obliczaj�ca wag� w�z�a,                               *)
    (* mvs - lista ruch�w do rozpatrzenia                             *)
    (* (beta, mb) - do tej pory obliczona waga i realizuj�cy j� ruch *)
    (* (waga sytuacji po * who, ruch).                               *)
    let choose mvs = 
      try 
	fold
	  (fun (beta, mb) m -> 
	    let (w, _) = alfabeta (move s who m) (other who) beta
            in 
	      if w * who >= alpha * who then raise (Obciecie (w, m)) else
	      if w * who > beta * who then (w, m) else (beta, mb))
	  (no_result who, pass)
	  mvs
      with Obciecie (beta, mb) -> (beta, mb)
    in 
      choose (moves s who);;
      
let game s = 
  let (_, m) = alfabeta s i (no_result (other i))
  in m;;


(***** Wersja ze spami�tywaniem. *****)
#use "spamietywanie.ml";;

let tab = ref empty_map;;

let rec alfabeta s who alpha = 
  memoize tab (function (s, who, alpha) -> 
    if finished s who then 
      (* Gra sko�czona *)
      (result s who, pass)
    else
      (* Gra trwa *)
      
      (* Funkcja obliczaj�ca wag� w�z�a,                               *)
      (* mvs - lista ruch�w do rozpatrzenia                            *)
      (* (beta, mb) - do tej pory obliczona waga i realizuj�cy j� ruch *)
      (* (waga sytuacji po * who, ruch).                               *)
      let choose mvs = 
	try 
	  fold
	    (fun (beta, mb) m -> 
	      let (w, _) = alfabeta (move s who m) (other who) beta
              in 
	        if w * who >= alpha * who then raise (Obciecie (w, m)) else
		if w * who > beta * who then (w, m) else (beta, mb))
	    (no_result who, pass)
	    mvs
	with Obciecie (beta, mb) -> (beta, mb)
      in 
        choose (moves s who))
    (s, who, alpha);;
      
let game s = 
  let (_, m) = alfabeta s i (no_result (other i))
  in m;;

