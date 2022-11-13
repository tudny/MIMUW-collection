
open List;;

let normalize lista =
  let rec aux li mini =
    match li with
    | [] -> []
    | (h :: t) ->
      let minimalny = min mini h in
      minimalny :: (aux t minimalny)
  in aux lista max_int

(* normalize [5; 6; 4; 3; 6; 2; 3];; *)

exception NoMiejsce

let rec odloz krazek pozostale =
  match pozostale with
  | []     -> raise NoMiejsce
  | h :: t ->
    if h >= krazek then
      t
    else
      odloz krazek t

(* odloz 3 [2; 2; 3; 3; 4; 5; 5];; *)

let krazki srednice krazki_li =
  let odw_srednice = rev (normalize srednice) in
  let rec aux sr_left kr_left =
    match kr_left with
    | []     -> length sr_left + 1
    | h :: t ->
      try aux (odloz h sr_left) t with
      | NoMiejsce -> 0
  in aux odw_srednice krazki_li;;

krazki [5; 6; 4; 3; 6; 2; 3] [3; 2; 5];;

krazki [3; 2; 1] [2; 4];;
