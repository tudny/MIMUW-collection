
(* Nie mam czasu, wiec robie na pale *)

let trojki l =
  let rec aux a b c acc =
    match a with
    | [] -> acc
    | x::ta ->
      match b with
      | [] -> aux ta ta ta acc
      | y::tb ->
        if a = b then aux a tb tb acc else
          match c with
          | [] -> aux a tb tb acc
          | z::tc ->
            if b = c then aux a b tc acc else
            if x + y <= z then
              aux a tb tb acc
            else
              aux a b tc ((x, y, z)::acc)
  in List.rev(aux l l l []);;

(* - pawlewicz *)

trojki [1; 2; 3; 4; 5; 6; 7; 8];;
