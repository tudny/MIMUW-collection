
type 'a tree = 
  | Node of 'a tree  * 'a * 'a tree 
  | Null
;;


(* exception Exit
let ultraleft t = 
  let rec loop t h acc = 
    match t with
    | Null -> if h < acc then raise Exit else h
    | Node (l, v, r) -> 
      loop l (h + 1) (loop r (h + 1) acc) in
  try let _ = loop t 0 0 in true 
    with
    | Exit -> false


let t = Node (Node (Null, 1, Null), 1, Node (Null, 1, Null));;
ultraleft t;;

let t = Node (t, 1, Null);;
ultraleft t;;

let t = Node (Null, 1, t);;
ultraleft t;;

ultraleft (Node ( Node ( Node (Null, 1, Null), 2, Null ), 4, Node (Node (Node (Node (Node (Null, 9, Null), 8, Null), 7, Null), 5, Null), 6, Null)));; *)



let ultraleft t = 
  let rec loop t h (wysokosc_ostatniego_liscia, czy) = 
    match t with
    | Null -> (h, czy && h >= wysokosc_ostatniego_liscia)
    | Node (l, v, r) -> 
      let (nowa_wysokosc_ostatnio_rozpatrywanego_liscia, nowe_czy) = loop r (h + 1) (wysokosc_ostatniego_liscia, czy) 
      in loop l (h + 1) (nowa_wysokosc_ostatnio_rozpatrywanego_liscia, nowe_czy) in
  snd (loop t 0 (0, true))

  let t = Node (Node (Null, 1, Null), 1, Node (Null, 1, Null));;
  ultraleft t;;
  
  let t = Node (t, 1, Null);;
  ultraleft t;;
  
  let t = Node (Null, 1, t);;
  ultraleft t;;
  
  ultraleft (Node ( Node ( Node (Null, 1, Null), 2, Null ), 4, Node (Node (Node (Node (Node (Null, 9, Null), 8, Null), 7, Null), 5, Null), 6, Null)));;
  
  