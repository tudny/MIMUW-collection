(* 
   Dane jest drzewo binarne, w postaci zgodnej z nast�puj�c� deklaracj�: 
*)
type drzewo =
    Ciekawe of (drzewo * drzewo) |
    Nudne of (drzewo * drzewo) |
    Puste;;

(* 
   Napisz procedur� �cie�ka: drzewo -> int, kt�ra dla 
   danego drzewa oblicza d�ugo�� najkr�tszej �cie�ki zaczynaj�cej si� 
   w korzeniu i przechodz�cej przez wszystkie ciekawe wierzcho�ki. 
*)


(* Rozwiazanie dynamiczne *)
let sciezka d = 
  (* Procedura pomocnicza. Dla ka�dego poddrzewa oblicza tr�jk�:
     (minimalna d�. �cie�ki schodz�cej z korzenia, 
     minimalna d�. �cie�ki schodz�cej i wracaj�cej do korzenia, 
     Czy w poddrzewie s� jakie� ciekawe wierzcho�ki). *)
  let rec pom drz = 
    match drz with
      Puste -> (0, 0, false) |
      _ ->
	let (l, p, ciekawe) = 
	  match drz with 
	    Ciekawe (l, p) -> (l, p, true) |
	    Nudne (l, p) -> (l, p, false)
	in
          let (l1, l2, lc) = pom l
	  and (p1, p2, pc) = pom p
	  in
	    if (not lc) && (not pc) then (0, 0, ciekawe) else
	    if not lc then (p1+1, p2+2, true) else
	    if not pc then (l1+1, l2+2, true) else
	    (min (l2 + p1 + 3) (l1 + p2 + 3), l2 + p2 + 4, true) 
  in
    let (x, _, _) = pom d
    in x;;
