
open Origami;;

let seed =
  try
    int_of_string Sys.argv.(1)
  with
  | Invalid_argument exp ->
    print_string "No argument. Setting the seed := 1201.\n";
    1201
  | Failure exp ->
    failwith "The argument should be castable to an int!"
;;


let act = ref seed;;
let rand () =
  let p = 7919 in
  let m = 1000000009 in
  act := !act * p mod m;
  !act
;;

let rand_int a =
  rand () mod a
;;

let rand_int_range a b =
  rand_int (b - a + 1) + a
;;

let cut_precision a =
  let precission = 1000000. in
  (float_of_int ( int_of_float (a *. precission) ) ) /. precission;;

let generate_cord () =
  let a = rand_int 10000000
  and b = rand_int_range 1 1000001
  in cut_precision ((float_of_int a) /. (float_of_int b))
;;

let generate_point () =
  (generate_cord (), generate_cord ())
;;

let generate_point_given (a, b) =
  let (c, d) = generate_point () in
  (a +. c, d +. b)
;;

let generate_line () =
  let p1 = generate_point ()
  and p2 = generate_point ()
  in (p1, p2)
;;

let generate_lines n =
  let rec gen n acc =
    if n = 0 then acc else
      let line = generate_line () in
      gen (n - 1) (line :: acc)
  in gen n []
;;

let string_of_point (a, b) =
  "(" ^ (string_of_float a) ^ ", " ^ (string_of_float b) ^ ")"
;;

let string_of_pair_of_points (p, q) =
  "(" ^ (string_of_point p) ^ ", " ^ (string_of_point q) ^ ")"
;;

let string_of_list_of_pairs_of_points li =
  let rec gen li acc =
    match li with
    | [] -> acc ^ "]"
    | h :: [] -> gen [] (acc ^ (string_of_pair_of_points h))
    | h :: t -> gen t (acc ^ (string_of_pair_of_points h) ^ "; ")
  in gen li "["
;;

type point = float * float;;

type start =
  | Prostokat of point * point
  | Kolko of point * float
;;

let save_to_file start lines name =
  let file = "tests/test" ^ name ^ ".ml" in
  let oc = open_out file in
  let fprint = Printf.fprintf oc in
  fprint "open Origami;;\n\n";
  fprint "let a = ";
  let pocz =
    match start with
    | Prostokat (p1, p2) -> Printf.fprintf oc "prostokat %s %s;;\n" (string_of_point p1) (string_of_point p2); prostokat p1 p2
    | Kolko (p, r) -> Printf.fprintf oc "kolko %s %s;;\n" (string_of_point p) (string_of_float r); kolko p r
  in
  Printf.fprintf oc "let lista = %s;;\n" (string_of_list_of_pairs_of_points lines);
  Printf.fprintf oc "let zlozenie = skladaj lista a;;\n";
  let zlozone = skladaj lines pocz in
  let n = (rand_int_range 1 100000) in
  for i = 0 to n do
    let p = generate_point () in
    let wynik = zlozone p in
    Printf.fprintf oc "assert ( (zlozenie %s) = (%d) );;\n" (string_of_point p) (wynik)
  done;
  close_out oc;
;;

let generate name =
  Printf.printf "Generowanie %s\n" name;
  let (lines, start) =
    let nr_of_lines = rand_int_range 0 2 in
    let lines = generate_lines nr_of_lines in
    let start =
      match rand_int 2 with
      | 0 -> let p = generate_point () in Prostokat (p, generate_point_given p)
      | _ -> Kolko (generate_point (), generate_cord ())
    in (lines, start)
  in
  save_to_file start lines name
;;

generate "Pierwszy";;
generate "Drugi";;
generate "Trzeci";;
generate "Czwarty";;
generate "Piaty";;

