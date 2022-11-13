
let tryw n = 
	let rec czy_zawiera k =
		if k >= n - 1 then false
		else if k * k mod n = 1 then true
		else czy_zawiera (k + 1)
	in czy_zawiera 2;;

tryw 2;;

