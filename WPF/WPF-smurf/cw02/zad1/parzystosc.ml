
let rec parzystosc x =
	if x mod 2 = 0 then 1 + parzystosc (x / 2)
	else 0;;
	
	
print_int (parzystosc 16);;
print_endline "";

print_int (parzystosc 5);;
print_endline "";

print_int (parzystosc 24);;
print_endline "";

