
let rec dub lista = 
	match lista with
	| [] 	 -> []
	| h :: t -> h :: h :: (dub t);;
	
dub [1; 2; 3; 4];;

