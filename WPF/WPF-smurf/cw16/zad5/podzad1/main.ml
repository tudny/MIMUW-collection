
type 'α tree =
  Node of 'α * 'α tree * 'α tree * 'α tree ref |
  Null;;

let infix t =
  let in_fix = ref [] in
  let rec walk trr =
    match trr with
    | Null -> ()
    | Node (v, l, r, p) ->
      walk l;
      in_fix := trr :: !in_fix;
      walk r;
  in
  walk t; !in_fix
;;

let fastyguj t =
  let in_fix = ref (infix t) in
  let rec walk trr =
    match trr with
    | Null -> ()
    | Node (v, l, r, p) ->
      walk l;
      match !in_fix with
      | [] -> invalid_arg "Coś poszło nie tak"
      | hd :: tl ->
        p := hd;
        in_fix := tl;
        walk r;
  in
  walk t
;;


