
 type 'a tree = Node of 'a  * 'a tree list;;


let rosnaca t = 
  let max (a, b) (c, d) = 
    if a > c then (a, b) else (c, d) in
  let rec aux (Node (v, lst)) = 
    List.fold_left
    (fun (maks_acc, lst_acc, maks_ja, lst_ja) (maks_spod_syna, lst_spod_syna, maks_z_synem, lst_z_synem) -> 
      match lst_z_synem with
      | [] -> invalid_arg "Miał być z synem!"
      | syn :: rest_syna ->
        let (ojcowa_dl, ojcowa_lst) = 
          if v < syn then 
            (maks_spod_syna + 1, v :: lst_z_synem)
          else
            (1, [v])
        in
          let (n_maks, n_lst) = max (maks_acc, lst_acc) (maks_spod_syna, lst_spod_syna) in
          let (n_maks, n_lst) = max (n_maks, n_lst) (ojcowa_dl, ojcowa_lst) in
          let (n_ojcowa, n_ojcowa_lst) = max (ojcowa_dl, ojcowa_lst) (maks_ja, lst_ja) in
          (n_maks, n_lst, n_ojcowa, n_ojcowa_lst)
    )
    (1, [v], 1, [v])
    (List.map (aux) lst) 
  in let (_, score, _, _) = aux t in score;;


let t = Node (5, [Node(4, [Node(5, [Node(6, [])])]); Node(6, []); Node(7, [])]);;
rosnaca t;;

let t = Node (5, [Node(4, [Node(5, [Node(6, [])])]); Node(6, [Node(10, [Node(20, [Node(70, [Node(10, [])]); Node(30, [Node(40, [])])])])]); Node(7, [])]);;
rosnaca t;;

let t1 = Node(30, [Node(40, [])]);;
rosnaca t1;;

let t2 = Node(20, [Node(70, [Node(10, [])]); Node(30, [Node(40, [])])]);;
rosnaca t2;;





