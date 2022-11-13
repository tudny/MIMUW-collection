(* Algorytm mini-max. *)

(* Gracze: max = 1, min = -1. *)
let i = 1;;

(* Przeciwnik gracza who *)
let other who = - who;;

(* Wyniki gier. *)
let won = 1;;
let lost = -1;;

(* Przyk³adowa gra. *)
#use "gra-w-kamienie.ml";;

open List;;

(* Funkcja wybieraj±ca z ci±gu par (w, s) parê o maksymalnym w. *)
let choose l = 
  match l with 
    [] -> failwith "Pusta lista ruchów!" |
    (h::t) -> 
      fold_left (fun (w, m) (w1, m1) -> if w1 > w then (w1, m1) else (w, m)) 
	   h t;;

(* Algorytm mini-max. *)
(* Ocena sytuacji s, gdy pierwszy ruch wykonuje who. *)
(* Wynikiem jest para (waga, ruch realizuj±cy wagê). *)
let rec minimax s who = 
  if finished s who then 
    (* Gra skoñczona *)
    (result s who, pass)
  else
    (* Gra trwa *)
    (* Funkcja przyporz±dkowuj±ca ruchowi parê *)
    (* (waga sytuacji po * who, ruch). *)
    let make_move m = 
      let (w, _) = minimax (move s who m) (other who)
      in (who * w, m)
    in 
      let
        (* waga optymalnego ruchu * who i optymalny ruch *)
        (w, m) = choose (map make_move (moves s who))
      in 
        (w * who, m);;

(* Optymalny ruch *)
let game s = 
  let (_, m) = minimax s i
  in m;;
