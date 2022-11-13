
let sumy li =
  let create (lst, acc) ele =
    let value = ele + acc in
    (value :: lst, value) in
  List.rev (fst (List.fold_left (create) ([], 0) li));;



sumy [1; 5; 2; 7; 12; 10; 5];;

assert (sumy [1; 5; 2; 7; 12; 10; 5] = [1; 6; 8; 15; 27; 37; 42])