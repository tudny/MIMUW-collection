
let czy_pierwsza n =
	if n < 2 then false else
	let rec check id n =
		if id = n then true
		else if n mod id = 0 then false
		else check (id + 1) n
	in check 2 n;;

assert ((czy_pierwsza 0) = false);;
assert ((czy_pierwsza 1) = false);;
assert ((czy_pierwsza 2) = true);;
assert ((czy_pierwsza 3) = true);;
assert ((czy_pierwsza 4) = false);;
assert ((czy_pierwsza 5) = true);;
assert ((czy_pierwsza 6) = false);;
assert ((czy_pierwsza 7) = true);;
assert ((czy_pierwsza 8) = false);;
assert ((czy_pierwsza 9) = false);;
assert ((czy_pierwsza 10) = false);;
assert ((czy_pierwsza 11) = true);;
assert ((czy_pierwsza 12) = false);;



