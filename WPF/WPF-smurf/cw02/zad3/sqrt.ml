
let sqrt n = 
	let rec aux x sum cnt =
		if sum > x then cnt - 1
		else aux x (sum + cnt + cnt - 1) (cnt + 1)
	in aux n 1 1;;
	
assert ((sqrt 1) = 1);;
assert ((sqrt 2) = 1);;
assert ((sqrt 3) = 1);;
assert ((sqrt 4) = 1);;
assert ((sqrt 5) = 1);;
assert ((sqrt 6) = 1);;
assert ((sqrt 7) = 1);;
assert ((sqrt 8) = 1);;
assert ((sqrt 9) = 1);;
assert ((sqrt 10) = 1);;

