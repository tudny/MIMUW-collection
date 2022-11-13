(* Obliczanie NWD algorytmem Euklidesa *)

(* Przez odejmowanie *)
let rec nwd x y =
  if x = y then x else
    if x > y then 
      nwd (x - y) y
    else
      nwd x (y - x);;

(* Przez dzielenie *)
let nwd x y = 
  let rec e x y = 
    if y = 0 then x else e y (x mod y)
  in
    if x > y then e x y else e y x;;

(* Przez parzysto¶æ *)
let parzyste x = x mod 2 = 0;;
let nwd x y =
  let rec pom x y a =  
    if x = y then a * x 
    else if parzyste x && parzyste y then pom (x / 2) (y / 2) (2 * a)
    else if parzyste x then pom (x / 2) y a
    else if parzyste y then pom x (y / 2) a
    else if x > y then pom (x - y) y a else pom x (y - x) a
  in
    pom x y 1;;
  


