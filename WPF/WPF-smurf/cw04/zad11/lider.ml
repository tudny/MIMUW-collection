
exception Empty_List;;

let lider lista =
  let rec aux li akt_lider cnt =
    match li with
    | []       -> akt_lider
    | (h :: t) ->
      if h = akt_lider then aux t akt_lider (cnt + 1)
      else if cnt = 0 then aux t h 1
      else aux t akt_lider (cnt - 1)
  in
  match lista with
    | []             -> raise Empty_List
    | (head :: tail) -> aux tail head 1;;
