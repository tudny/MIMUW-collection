(* Logarytm całkowitoliczbowy. *)

let square x = x * x;;

(* T(n) = O(log n),  M(n) = O(log n) *)
let rec int_log n = 
  if n = 1 then 0 
  else 1 + int_log (n / 2);;

(* T(n) = O(log n),  M(n) = O(1) *)
let int_log n = 
  let rec pom n a = 
    if n = 1 then a 
    else pom (n / 2) (a + 1)
  in pom n 0;;

(* T(n) = O((log log n)^2),  M(n) = O(1) *)
let int_log n = 
  let rec pom i j k m acc = 
    (* j = 2^i, k = 2^j, k <= m, wynik = int_log m + acc = int_log n *)
    if m = 1 then acc
    else if square k > m then  
	pom 0 1 2 (m / k) (acc + j)
    else 
      pom (i+1) (2*j) (square k) m acc
  in pom 0 1 2 n 0;;

(* T(n) = O(log log n),  M(n) = O(log log n) *)
let int_log n = 
  let rec gen i j k acc = 
    (* j = 2^i, k = 2^j, sqrt k <= n, acc = [(2^(i-1), 2^(2^(i-1))); \dots (2, 2^2); (1, 2^1)] *)
    if k > n then acc
    else gen (i+1) (2*j) (square k) ((j, k) :: acc)
  and scan l m = 
    (* l = [(2^i, 2^(2^i)); \dots (2, 2^2); (1, 2^1)],  m < 2^(2^(i+1))) *)
    match l with 
      []        -> 0 |
      (j, k)::t -> 
	if k <= m then 
	  j + scan t (m / k)
	else
	  scan t m
  in
    scan (gen 0 1 2 []) n;;

(* T(n) = O(log log n),  M(n) = O(1) *)
let int_log n = 
  let rec gen i j k l m = 
    (* j = Fib_i, k = Fib_(i+1), l = 2^Fib_i, m = 2^Fib_(i+1),  l <= n, *)
    if m > n then 
      scan i j k l m n 0
    else 
      gen (i+1) k (j+k) m (l*m)
  and scan i j k l m nn acc = 
    (* j = Fib_i, k = Fib_(i+1), l = 2^Fib_i, m = 2^Fib_(i+1),  m > nn, int-log n = acc + int-log nn *)
    if nn = 1 then acc 
    else if l <= nn then 
      scan (i-1) (k-j) j (m/l) l (nn/l) (acc+j)
    else
      scan (i-1) (k-j) j (m/l) l nn acc
  in
    gen 1 1 1 2 2;;


assert (int_log 1 = 0);;
assert (int_log 42 = 5);;
assert (int_log 16 = 4);;
assert (int_log 256 = 8);;
assert (int_log 1024 = 10);;
