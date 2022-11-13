
let codrugi li =
  let brac_nie_brac (lst, boo) ele =
    if boo then (ele :: lst, false)
    else (lst, true) in
  List.rev (fst (List.fold_left (brac_nie_brac) ([], false) li));;


assert (codrugi [1; 2; 3; 4; 5] = [2; 4]);;
assert (codrugi [] = []);;
assert (codrugi [1] = []);;
assert (codrugi [1; 2] = [2]);;