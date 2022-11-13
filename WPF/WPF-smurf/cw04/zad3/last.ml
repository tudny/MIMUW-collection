
exception EmptyList of string;;

let rec last lista = 
	match lista with
	| [] 	  -> raise (EmptyList "The list is empty!")
	| h :: [] -> h
	| h :: t  -> last t;;
	
	
	
last [1; 2; 3; 4];;

