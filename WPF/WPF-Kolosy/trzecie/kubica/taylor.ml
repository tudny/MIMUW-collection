(* Aproksymacja funkcji za pomoc± szeregów Taylora. *)


(* Sk³adanie funkcji itp. *)
let id x = x;;

let compose f g x = f (g x);;

let rec iterate n f =
  if n = 0 then id else compose (iterate (n-1) f) f;;

(* Pochodna *)
let epsilon = 0.0001;;

let rozniczka f x dx = (f (x +. dx) -. f x) /. dx;;

let pochodna f x = rozniczka f x epsilon;;


(* Suma czê¶ciowa szeregu. *)
let rec szereg f n = 
  if n = 0 then 0.0 else f n +. szereg f (n-1);;


(* Silnia *)
let rec silnia n = 
  if n < 2 then 1 else n * silnia (n - 1);;


(* Rozwi±zanie nieoptymalne, ale proste. *)
let potega x n = iterate n (fun a -> a *. x) 1.0;;

let wyraz f x0 x n = 
  iterate (n-1) pochodna f x0 /. float (silnia (n-1)) *. potega (x -. x0) (n - 1);;

let taylor f x0 n x = szereg (wyraz f x0 x) n;;

(* Testy *)
let w x = 3. *. x *. x *. x -. 2. *. x *. x +. 1.;;
let v = taylor w 0. 4;;

List.map w [0.;1.;2.;3.];;
(* = [1.; 2.; 17.; 64.] *)

List.map v [0.;1.;2.;3.];;
(* = [1.; 2.00074467938804368; 17.003557286119463; 64.0087057630468479] *)
