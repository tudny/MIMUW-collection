
let czy_pal lista = lista = List.rev lista;;

let palindrom li = 
let rec in_fiks lista acc = 
  let rec check li acc = 
    match li with
    | [] -> acc
    | h :: t -> if czy_pal li 
      then max acc (List.length li)
      else check t acc
  in match lista with
  | [] -> acc
  | h :: t -> in_fiks t (check (List.rev lista) acc) 
in in_fiks li 0;;

palindrom [10; 20; 30; 20; 1; 4; 2; 3; 3; 2; 4; 5];;
