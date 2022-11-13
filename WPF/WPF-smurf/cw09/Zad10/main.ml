
let widoczne lista =
  let suma = List.fold_left (+) 0 lista in
  let (w, _) =
    List.fold_left
      (fun (acc, s_p) ele ->
         let s_z = suma - s_p - ele in
         if s_z = s_p then
           (ele :: acc, s_p + ele)
         else
           (acc, s_p + ele)
      )
      ([], 0) lista
  in List.rev w;;

widoczne [2; 1; -1; 1; -1];;
