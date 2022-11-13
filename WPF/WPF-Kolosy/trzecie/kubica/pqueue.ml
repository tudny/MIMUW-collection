type 'a pri_queue =
 Node of 'a * 'a pri_queue * 'a pri_queue * int | Leaf
exception EmptyQueue
let empty_queue = Leaf
let is_empty q = q = Leaf

let size q =
 match q with
 Leaf -> 0 |
 Node (_, _, _, n) -> n

let get_root h =
 match h with
 Leaf -> raise EmptyQueue |
 Node (r, _, _, _) -> r

let set_root h r =
 match h with
 Leaf -> Node (r, Leaf, Leaf, 1) |
 Node (_, l, p, n) -> Node (r, l, p, n)

let rec put h x =
(*
 let rotate r l p =
  let n = size l + size p + 1
  in match (l,p) with
   (Leaf, Leaf) -> Node (r, Leaf, Leaf, n) |
   (_, Leaf) -> if r >= get_root l then Node (r, l, Leaf, n)
                else Node (get_root l, set_root l r, Leaf, n)|
   (Leaf, _) -> if r >= get_root p then Node (r, Leaf, p, n)
                else Node (get_root p, Leaf, set_root p r, n) |
   _ -> let rl = get_root l
        and rp = get_root p
        in let (r1, rl1, rp1) =
         if r >= rl then
          if r >= rp then (r, rl, rp) else (rp, rl, r)
         else
          if rl >= rp then (rl, r, rp) else (rp, rl, r)
         in Node (r1, set_root l rl1, set_root p rp1, n)
 in
(* Niepotrzebne *)
*)
 match h with
 Leaf -> Node (x, Leaf, Leaf, 1) |
 Node (r, l, p, n) ->
  if(size l <= size p)
   then Node((max x r), (put l (min x r)), p, (n+1)) (* Mozna tak. *)
   else Node((max x r), l, (put p (min x r)), (n+1))

let getmax h =
 let rec del h =
  match h with
  Leaf -> h |
  (* Node (_, Leaf, Leaf, _) -> Leaf |  (* niepotrzebne *) *)
  Node (_, l, Leaf, _) -> l |
  Node (_, Leaf, p, _) -> p |
  Node (_, (Node (rl, _, _, _) as l), (Node (rp, _, _, _) as p), n) ->
   if rl >= rp then
    Node (rl, del l, p, n - 1)
   else
    Node (rp, l, del p, n - 1)
 in (get_root h, del h)
