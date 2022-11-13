(* Mnożenie metodą rosyjskich chłopów *)

let parzyste x = x mod 2 = 0;;

let rec razy x y = 
  let rec pom x1 y1 a =
    assert ((x1 >= 0) && (x1 * y1 + a = x * y));
    if x1 = 0 then 
      a 
    else if parzyste x1 then 
      pom (x1 / 2) (2 * y1) a
    else
      pom (x1 - 1) y1 (a + y1)
in
  if abs x > abs y then 
    razy y x
  else if x > 0 then 
    pom x y 0 
  else
    pom (-x) (-y) 0;;

assert (razy 2 2 = 4);;
assert (razy 5 (-3) = -15);;
assert (razy 42 0 = 0);;
