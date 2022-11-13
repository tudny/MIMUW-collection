
(* Rotacja z dodatkową pamięcią *)

let rotacja arr k =
  let n = Array.length arr in
  let dup = Array.init n (fun i -> arr.(i)) in
  Array.iteri (fun i ele -> arr.(i) <- dup.((i - k + n) mod n)) arr
;;

let arr = [|1; 2; 3; 4; 5; 6; 7|];;

rotacja arr 3;;
arr;;

(* Da się bez na odwracanie tablicy *)
(* Rot(a, b) -> rotacja tablicy na przedziale (a, b) *)
(* Rot(1, n) *)
(* Rot(1, k) *)
(* Rot(k + 1, n) *)


