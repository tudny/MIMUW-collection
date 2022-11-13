(* Gra w kamienie - czê¶æ zale¿na od gry *)
(* Stan gry = liczba pozosta³ych kamieni *)
(* Ruch = liczba zdejmowanych kamieni. *)

#use "listy.ml";;

type board = int;;
type move = int;;

(* Funkcja generuj±ca dla danej sytuacji s i gracza who listê ruchów. *)
let moves s who = 
  match s with 
    0 -> [] |
    1 -> [1] |
    2 -> [1; 2] |
    _ -> [1; 2; 3];;
    
(* Pusta lista ruchów, to koniec gry. *)
let finished s who = moves s who = [];;

(* Pas - brak ruchu, np. gdy gra skoñczona. *)
let pass = 0;;

(* Funkcja przypisuj±ca wagi li¶ciom, s = sytuacja, who = gracz. *)
(* Zak³adamy, ¿e finished s who. *)
let result s who = 
  if finished s who then 
    if who = i then lost else won
  else failwith "The game is not finished yet!";;

(* Rezultat wykonania ruchu m w sytuacji s przez gracza who *)
let move s who m = s - m;;

