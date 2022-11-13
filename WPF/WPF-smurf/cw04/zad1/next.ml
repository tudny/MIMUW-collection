
let next x = 
	let rec aux n acc = 
		if n = 0 then acc
		else aux (n - 1) (n :: acc)
	in aux x [];;
	
next 10;;

