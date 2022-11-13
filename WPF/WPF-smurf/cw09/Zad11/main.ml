
let podzial lista =
  if lista = [] then [] else
    let pier = List.hd (List.rev lista) in
    let (_, acc, big_acc) =
    List.fold_right
      (fun ele (poprz, acc, big_acc) ->
         if ele = poprz then (ele, [ele], acc :: big_acc)
         else (ele, ele :: acc, big_acc)
      )
      (List.rev (List.tl (List.rev lista)))
      (pier, [pier], [])
    in if acc = [] then big_acc else acc :: big_acc;;

podzial [3; 2; 2; 5; 7; 5; 4; 4; 3; 1];;
