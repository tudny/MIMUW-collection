
let compose f g = fun x -> f (g x)

let rec iterate n f = 
  if n = 0 then
    (fun x -> x)
  else 
    compose f (iterate (n - 1) f);;


iterate 2 (function x -> x * (x + 1)) 2;;
iterate 3 (function x -> x * (x + 1)) 1;;

