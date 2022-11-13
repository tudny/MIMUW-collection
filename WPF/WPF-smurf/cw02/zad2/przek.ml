

let przeksztalc n = 
	let rec przek n acc = 
		if n = 0 then acc
		else przek (n / 10) ((n mod 10) + (acc * 10))
	in przek n 0;;
	
	
print_int(przeksztalc 1234);;
	
