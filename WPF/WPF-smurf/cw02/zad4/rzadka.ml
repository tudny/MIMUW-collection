
let rec czy_rzadka x =
	if x < 3 then true
	else if x mod 4 = 3 then false
	else czy_rzadka (x / 2);;

let rec rzadkie n =
	if czy_rzadka (n + 1) then n + 1
	else rzadkie (n + 1);;



print_int (rzadkie 42);;

