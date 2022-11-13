
let append li ele =
  List.fold_right (fun el acc -> el :: acc) li [ele]

let lista = [1; 2; 3; 4];;
append lista 5;;
