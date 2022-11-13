
let rec convert l =
	let transform (n, x) = 
		let rec loop id = 
			if id = 0 then []
			else x :: (loop (id - 1))
		in loop n
	in
	match l with
	| [] 	 -> []
	| h :: t -> (transform h) @ (convert t);;
	
convert [(1, 2); (2, 3); (4, 10)];;
