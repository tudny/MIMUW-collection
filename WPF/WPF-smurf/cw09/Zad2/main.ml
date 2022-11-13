
let neg pre x = not (pre x)

let exists p li =
  List.length (List.filter (p) li) > 0


(* Prawo De Morgana *)
(* not forall p = exists neg p *)
(* forall p = not exists neg p *)

let forall pred li =
  not (exists (neg pred) li)


let lista = [1; 2; 3; 4; 5; 6; 7; 8];;
let predykat = fun x -> x < 9;;

forall (predykat) lista;;

