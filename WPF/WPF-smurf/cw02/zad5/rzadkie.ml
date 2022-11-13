
let int_of_bool b = 
	if b then 1 else 0;;


let rec czy_rzadka x = 
	if x < 3 then true
	else if x mod 4 = 3 then false
	else czy_rzadka (x / 2);;
		

let rzadkie n = 
	let rec aux id up_to acc =
		 if id - 1 = up_to then acc
		 else aux (id + 1) up_to (acc + int_of_bool (czy_rzadka id))
	in aux 1 n 0;;
	

rzadkie 42;;

