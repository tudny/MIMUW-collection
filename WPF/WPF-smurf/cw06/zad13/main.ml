
type 'a drzewo = Node of 'a * 'a drzewo * 'a drzewo | Null

type 'a drzewoSized = NodeS of 'a * 'a drzewoSized * 'a drzewoSized * int | NullS

let polowa t = 
  let size t = 
    match t with
    | NullS -> 1
    | NodeS (_, _, _, s) -> s
  in
  let rec convert t = 
    match t with 
    | Null -> NullS
    | Node (v, l, r) ->
      let l = convert l 
      and r = convert r in 
      let s_l = size l
      and s_r = size r in
      NodeS (v, l, r, s_l + s_r)
    in
  let t = convert t in
  let superS = size t in
  let rec loop acc t =
    match t with 
    | NullS -> acc
    | NodeS (_, l, r, s) ->
      let acc = loop acc l in 
      let acc = loop acc r in
      if 2 * s >= superS && s < size acc then t else acc
  in loop t t
;;




