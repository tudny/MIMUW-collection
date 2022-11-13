
let head l n =
	if n >= List.length l then l else
	let rec head_aux (h :: t) n =
		if n = 0 then []
		else h :: (head_aux t (n - 1))
	in head_aux l n;;
	
head [1; 2; 3; 4; 5; 6] 7;;



let tail l n = 
	head (List.rev l) n;;

tail [1; 2; 3; 4; 5; 6] 4;;

