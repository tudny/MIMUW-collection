(* Gra w kamienie - cz�� zale�na od gry *)
(* Stan gry = liczba pozosta�ych kamieni *)
(* Ruch = liczba zdejmowanych kamieni. *)

#use "listy.ml";;

type board = int;;
type move = int;;

(* Funkcja generuj�ca dla danej sytuacji s i gracza who list� ruch�w. *)
let moves s who = 
  match s with 
    0 -> [] |
    1 -> [1] |
    2 -> [1; 2] |
    _ -> [1; 2; 3];;
    
(* Pusta lista ruch�w, to koniec gry. *)
let finished s who = moves s who = [];;

(* Pas - brak ruchu, np. gdy gra sko�czona. *)
let pass = 0;;

(* Funkcja przypisuj�ca wagi li�ciom, s = sytuacja, who = gracz. *)
(* Zak�adamy, �e finished s who. *)
let result s who = 
  if finished s who then 
    if who = i then lost else won
  else failwith "The game is not finished yet!";;

(* Rezultat wykonania ruchu m w sytuacji s przez gracza who *)
let move s who m = s - m;;

