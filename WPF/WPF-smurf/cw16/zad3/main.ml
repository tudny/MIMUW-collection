
let duplikat arr =

  let go x = arr.(x) in

  let p1 = ref (go 0)
  and p2 = ref (go (go 0)) in

  while !p1 <> !p2 do
    p1 := go !p1;
    p2 := go (go !p2);
  done;

  let p1 = ref (go !p1)
  and p2 = ref (go (go !p1)) in

  let len = ref 1 in
  let inc a =
    a := !a + 1 in

  while !p1 <> !p2 do
    inc len;
    p1 := go !p1;
    p2 := go (go !p2);
  done;

  let rec go_n n p =
    if n = 0 then p
    else go_n (n - 1) (go p) in

  let p1 = ref (go_n !len 0)
  and p2 = ref 0 in

  while !p1 <> !p2 do
    inc len;
    p1 := go !p1;
    p2 := go !p2;
  done;

  !p1
;;

let arr = [|1; 2; 3; 2; 5; 4|];;
let cykl = duplikat arr;;
assert (cykl = 2);;

let arr = [|1; 2; 4; 5; 1; 3|];;
let cykl = duplikat arr;;
assert (cykl = 1);;

let arr = [|1; 2; 3; 4; 5; 5|];;
let cykl = duplikat arr;;
assert (cykl = 5);;

let arr = [|2; 1; 1; 1; 1; 1; 1; 1; 1|];;
let cykl = duplikat arr;;
assert (cykl = 1);;
