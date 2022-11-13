
type 'a tree = 
  | Node of 'a tree * 'a * 'a tree 
  | Nil
;;

let liscie t = 
  let rec loop t acc = 
    match t with
    | Nil -> acc
    | Node (Nil, v, Nil) -> v :: acc
    | Node (l, _, r) -> loop r (loop l acc)
  in loop t []
;;

let cnt = ref 0;;
let inc () = 
  let temp = !cnt in 
  cnt := !cnt + 1; temp
;;

let t () = Node (Nil, inc (), Nil);; 
let t1 () = Node (t (), inc (), t ());;

let t1 = Node (t1 (), inc (), t1 ());;

liscie t1;;