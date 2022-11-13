
let zlozenie lista =
  List.fold_right (fun z f -> (fun x -> z (f x))) lista (fun x -> x);;

let li = [(fun x -> 2 * x); (fun x -> x * x); (fun x -> x / 2)];;

let z = zlozenie li;;

z 10;;
(* - : int = 50 *)
(* z(x) = 2*((x/2)^2) *)


type ulamek = Zwykly of int * int |
              Dziesietny of float;;


let pi_zwykle = Zwykly (22, 7);;
let pi_dziesietne = Dziesietny (3.14);;

let get_float_value a =
  match a with
  | Zwykly (a, b) -> (float_of_int a) /. (float_of_int b)
  | Dziesietny a -> a;;

get_float_value pi_zwykle;;
get_float_value pi_dziesietne;;

get_float_value (Zwykly (21, 37));;
get_float_value (Dziesietny (4.20));;


let get_float_value = function
  | Zwykly (a, b) -> (float_of_int a) /. (float_of_int b)
  | Dziesietny a -> a;;

get_float_value pi_zwykle;;
get_float_value pi_dziesietne;;

get_float_value (Zwykly (21, 37));;
get_float_value (Dziesietny (4.20));;
