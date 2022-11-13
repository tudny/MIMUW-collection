open Origami;;

let centr = (0., 0.);;

let gen n =
  let rec aux acc i =
    if i > n then acc
    else
      aux (((float_of_int(1 lsl i), 0.),
            (float_of_int(1 lsl i), 1.))::acc) (i + 1)
  in
  aux [] 0;;

let const = 24;;

let l = gen const;;

let a = kolko centr (max_float /. 4.);;

let a = skladaj l a;;

assert(a centr = (1 lsl const) + 1);;
assert(a (1., 1.) = 1 lsl const);;
assert(a (-1., -1.) = 2);;
assert(a ((-.max_float) /. 2., (-.max_float) /. 2.) = 0);;




let x = prostokat (-16., -16.) (16., 16.);;
let a = (0., -16.);;
let b = (0., 16.);;
let c = (-16., 0.);;
let d = (16., 0.);;
let x = skladaj [(a,d);(d,b);(b,c);(c,a)] x;;

assert (x (0., 0.) = 5);;
assert (x (6., 0.) = 3);;
assert (x a = 1);;

