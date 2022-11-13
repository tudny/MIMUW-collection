

type 'a tree = 
  | Node of 'a tree * 'a * 'a tree
  | Null
;;

let potnij p t = 
  let rec loop acc = function
  | Null -> (acc, false)
  | Node (l, v, r) -> 
    let (acc_l, czy_l) = loop acc l in
    let (acc_r, czy_r) = loop acc_l r in 
    let nl = if czy_l then Null else l
    and nr = if czy_r then Null else r in
    let akt = Node (nl, v, nr) in
    if p v then (akt :: acc, true)
    else (acc, false)
  in fst (loop [] t)
;;

let pred a = abs a < 5
;;

let t1 = Node (Null, 1, Null);;
let t2 = Node (Null, 2, Null);;

potnij pred t1;;

let t1 = Node (t1, 10, t2);;



