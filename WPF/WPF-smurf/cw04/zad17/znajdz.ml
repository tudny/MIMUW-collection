


let znajdz matrix wartosc =
  let rec aux_k y x wiersz acc =
    match wiersz with
    | [] -> acc
    | h :: t ->
      let nowy_acc = if h = wartosc then ((x, y) :: acc) else acc
      in aux_k y (x + 1) t nowy_acc 
  in let rec aux_w y wiersze acc =
    match wiersze with
    | [] -> acc
    | h :: t -> aux_w (y + 1) t (aux_k y 1 h acc)
  in aux_w 1 matrix [];;

znajdz [[4; 3; 2; 1]; [8; 7; 3; 2]; [15; 10; 4; 1]; [40; 30; 20; 10]] 10;;
