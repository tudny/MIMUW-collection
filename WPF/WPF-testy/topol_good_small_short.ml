open Topol;;

exception False_exp
;;

exception True_exp
;;

let are_duplicates li = 
  let li = List.sort compare li in
  match li with
  | [] | [_] -> false
  | hd :: tl -> try let _ = List.fold_left
    ( fun last ele ->
      if ele = last then raise True_exp
      else ele
    ) hd tl in false
    with
    | True_exp -> true
;;

let remove_duplicates = function 
  | [] -> []
  | [x] -> [x]
  | hd :: tl ->
    List.fold_left ( fun (last, acc) ele ->
      if ele = last then (ele, acc)
      else (ele, ele :: acc)
    ) (hd, [hd]) tl |> snd |> List.rev
;;

let does_match data output = 
  let (first, second) = List.split data in
  let data = first @ List.flatten second in
  let data = List.sort compare data |> remove_duplicates
  and output = List.sort compare output in
  data = output
;;

let verify data output = 
  if are_duplicates output then false else
  if does_match data output |> not then false else

  let (map, _) = List.fold_left ( fun (_acc, iter) ele ->
      (PMap.add ele iter _acc, iter + 1)
    ) (PMap.empty, 0) output in

  let get_id ele = PMap.find ele map in

  let single (x, row) = 
    let x_id = get_id x in
    List.iter ( fun ele ->
      let ele_id = get_id ele in
      if x_id >= ele_id then raise False_exp
    ) row in 

  try
    List.iter ( fun ele -> 
      single ele
    ) data; true
  with
  | False_exp -> false
;;
assert ((let data = [(8, [4; 5; 7; 7; 9; 0]); (2, [4; 9; 5; 9; 6])] in let toped = Topol.topol data in verify data toped) );
assert ((let data = [(8, [2]); (0, [8; 8; 4])] in let toped = Topol.topol data in verify data toped) );
assert ((let data = [(3, [6; 9; 8; 7; 7; 0; 8; 8]); (3, []); (7, []); (2, [8; 0])] in let toped = Topol.topol data in verify data toped) );
assert ((let data = [(0, [4; 2]); (7, [4; 2; 5; 0])] in let toped = Topol.topol data in verify data toped) );
assert ((let data = [(7, []); (6, [0; 8; 8; 0; 9; 3; 4]); (5, [4; 8; 2])] in let toped = Topol.topol data in verify data toped) );
assert ((let data = [(5, [8; 3; 3; 9; 7; 0; 9; 0; 6]); (7, []); (6, [3]); (9, [2]); (4, [3; 7; 6])] in let toped = Topol.topol data in verify data toped) );
assert ((let data = [(5, [6; 1; 1; 9; 0]); (7, [2; 8; 4; 4; 6; 9; 8; 5]); (7, [8; 1; 9; 8; 8]); (4, [])] in let toped = Topol.topol data in verify data toped) );
assert ((let data = [(2, [4])] in let toped = Topol.topol data in verify data toped) );
assert ((let data = [(7, [1; 0; 9; 6; 0])] in let toped = Topol.topol data in verify data toped) );
assert ((let data = [(9, [3; 2; 3; 3; 2])] in let toped = Topol.topol data in verify data toped) );
assert ((let data = [(9, [2; 8; 1; 1; 2]); (3, [7; 6; 9; 5; 5; 2; 2]); (8, []); (5, [4; 9; 0]); (8, [])] in let toped = Topol.topol data in verify data toped) );
assert ((let data = [(2, [9; 0]); (0, [9])] in let toped = Topol.topol data in verify data toped) );
assert ((let data = [(6, [7; 7; 3; 3]); (1, []); (1, [4; 2; 9; 7; 4; 3; 3]); (5, [9; 4; 3; 7; 9; 0; 2; 8; 8])] in let toped = Topol.topol data in verify data toped) );
assert ((let data = [(8, []); (2, [])] in let toped = Topol.topol data in verify data toped) );
assert ((let data = [(4, [7; 1; 0; 3; 8; 6; 2])] in let toped = Topol.topol data in verify data toped) );
assert ((let data = [(2, [9; 0; 9; 4; 1]); (2, [0; 3; 6; 0; 6])] in let toped = Topol.topol data in verify data toped) );
assert ((let data = [(3, [6; 7; 7; 1; 8; 9; 4; 4])] in let toped = Topol.topol data in verify data toped) );
assert ((let data = [(2, []); (1, []); (5, [2]); (0, [5; 3; 3]); (4, [7; 9; 8; 1])] in let toped = Topol.topol data in verify data toped) );
