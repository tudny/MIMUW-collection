
let rec parzystosc x =
	if x mod 2 = 0 then 1 + parzystosc (x / 2)
	else 0;;

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

(* [(n_1, x_1); (n_2, x_2)] *)


let decomp x = 
  let parz = parzystosc x in
  let pot_dwa = 1 lsl parz in
  let rest = x / pot_dwa in
  (parz + 1, (rest + 1) / 2);;

let dekompresuj lista = 
  let rec aux li acc = 
    match li with
    | []   -> acc
    | h::t -> aux t ((decomp h) :: acc)
  in convert (List.rev (aux lista []));;

dekompresuj [1; 3; 3; 9; 42; 3];;
