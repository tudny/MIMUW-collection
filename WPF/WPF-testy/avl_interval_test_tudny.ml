
open ISet


(* Prosty test poprawności dodawania przedziałów *)
let a = empty;;
let a = add (-10, 10) a;;
let a = add (20, 30) a;;
let a = add (-100, -20) a;;

let elts = elements a;;
let elts_good = [(-100, -20); (-10, 10); (20, 30)];;
assert (elts = elts_good);;

let a = add (-200, 200) a;;
let elts = elements a;;
let elts_good = [(-200, 200)];;
assert (elts = elts_good);;

let a = add (min_int, min_int + 1) a;;
let a = add (max_int - 1, max_int) a;;
let elts = elements a;;
let elts_good = [(min_int, min_int + 1); (-200, 200); (max_int - 1, max_int)];;
assert (elts = elts_good);;


(* Generator par singletonów *)
let generate_singletons a b =
  let rec loop id acc =
    if id = b + 1 then acc
    else loop (id + 1) ((id, id) :: acc)
  in loop a []

let add_all = List.fold_left (fun a s -> add s a);;

let singletons = generate_singletons 1 20;;
let b = add_all empty singletons;;
let singletons = generate_singletons 40 100;;
let c = add_all b singletons;;
let singletons = generate_singletons min_int (min_int + 10);;
let d = add_all c singletons;;
let singletons = generate_singletons (max_int - 10) max_int;;
let e = add_all d singletons;;

let elts = elements e;;
let elts_good = [(min_int, min_int + 10); (1, 20); (40, 100); (max_int - 10, max_int)];;
assert (elts = elts_good);;


(* Generowanie przedziałów *)
let generate_intervals a b step size =
  let rec loop id acc =
    if id > b then acc
    else loop (id + step) ((id, id + size - 1) :: acc)
  in List.rev (loop a []);;


let intervals = generate_intervals 0 100 10 4;;
let f = add_all empty intervals;;

let elts = elements f;;
let elts_good = intervals;;
assert (elts = elts_good);;


let intervals = generate_intervals 4 104 10 6;;
let g = add_all f intervals;;

let elts = elements g;;
let elts_good = [(0, 109)];;
assert (elts = elts_good);;




let intervals = [(0, 0); (2, 3); (5, 7); (9, 12); (14, 18)];;
let h = add_all empty intervals;;

let elts = elements h;;
let elts_good = intervals;;
assert (elts = elts_good);;


let intervals = [(1, 1); (6, 10)];;
let i = add_all h intervals;;

let elts = elements i;;
let elts_good = [(0, 3); (5, 12); (14, 18)];;
assert (elts = elts_good);;




let interval = (6, 10);;
let j = remove interval h;;

let elts = elements j;;
let elts_good = [(0, 0); (2, 3); (5, 5); (11, 12); (14, 18)];;
assert (elts = elts_good);;



let interval = (3, 11);;
let k = remove interval j;;

let elts = elements k;;
let elts_good = [(0, 0); (2, 2); (12, 12); (14, 18)];;
assert (elts = elts_good);;
assert (mem 0 k);;
assert (mem 2 k);;
assert (mem 12 k);;
assert (mem 14 k);;
assert (mem 15 k);;
assert (mem 16 k);;
assert (mem 17 k);;
assert (mem 18 k);;
assert (not (mem 20 k));;
assert (not (mem 21 k));;
assert (not (mem (-1) k));;
assert (not (mem 1 k));;
assert (not (mem 13 k));;
assert (not (mem max_int k));;
assert (not (mem min_int k));;
assert ((below 10 k) = 2);;
assert ((below 100 k) = 8);;
assert ((below min_int k) = 0);;
assert ((below max_int k) = 8);;
assert ((below 15 k) = 5);;




let interval = (6, 100);;
let l = remove interval k;;

let elts = elements l;;
let elts_good = [(0, 0); (2, 2)];;
assert (elts = elts_good);;



let interval = (min_int, max_int);;
let m = remove interval l;;

let elts = elements m;;
let elts_good = [];;
assert (elts = elts_good);;
assert (is_empty m);;



let full = add (min_int, max_int) empty;;
assert (mem max_int full);;
assert (mem min_int full);;



let t1 = add (min_int, min_int + 100) empty;;
let t2 = add (max_int - 100, max_int) t1;;
assert ((below 0 t2) = 101);;
assert ((below max_int t2) = 202);;

let t3 = add (0, 0) t2;;
assert ((below 10 t3) = 102);;

