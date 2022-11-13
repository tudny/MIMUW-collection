
(* Punkt x, y *)
type point = int * int;;

(* Przedział *)
type range = int * int;;

(* Dwa skrajne punkty prostokąta *)
type rectangle = point * point;;


let axis ((n, m) : rectangle) which =
  ((which n, which m) : range);;

let intersection ((a, b) : range) ((c, d) : range) =
  ((max a c, min b d) : range);;

let merge ((a, b) : range) ((c, d) : range) =
  (((a, c), (b, d)) : rectangle);;


let przeciecie (a : rectangle) (b : rectangle) = 	merge (intersection (axis a fst) (axis b fst)) (intersection (axis a snd) (axis b snd));;



(* Testy *)

assert (przeciecie ((1, 3), (5, 6)) ((2, 1), (7, 5)) = ((2, 3), (5, 5)));;
assert (przeciecie ((-1, 1), (5, 2)) ((1, -1), (3, 4)) = ((1, 1), (3, 2)));;
let ((x1, y1), (x2, y2)) = przeciecie ((0, 0), (1, 1)) ((2, 1), (4, 4))
in assert (x1 > x2 || y1 > y2);;
