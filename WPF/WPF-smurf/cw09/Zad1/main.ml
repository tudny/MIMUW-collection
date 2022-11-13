
let exists p li =
  List.length (List.filter (p) li) > 0



let lista = [1; 2; 3; 4; 6; 11; 21];;
let predykat x = x mod 5 = 0;;

assert (exists (predykat) lista = false);;
