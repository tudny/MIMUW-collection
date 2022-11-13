
open Topol;;

let print_list li =
  List.iter (Printf.printf "%d, ") li; print_string "\n"
;;

let g = [(1, [2; 3]); (2, [4]); (3, [4; 5]); (4, []); (5, [])];;

let t = Topol.topol g;;

print_list t;;

let g2 = [(1, [2]); (2, [1])];;

assert (try let _ = topol g2 in false with _ -> true);;

let b = [(1, [2]); (2, [3]); (3, [4])];;

let b2 = Topol.topol b;;

print_list b2;;
