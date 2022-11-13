
exception True_exp
;;

let kwadrat (arr : bool array array) =
  let n = Array.length arr
  and m = Array.length arr.(0)
  in 
  let sum_pref =
    let pref = Array.make_matrix n m 0 in 
    let get x y = try pref.(x).(y) with _ -> 0 in
    Array.iteri ( fun x row -> 
        Array.iteri ( fun y el -> 
            pref.(x).(y) <- (get x (y - 1)) + (get (x - 1) y)
              - (get (x - 1) (y - 1)) + if not el then 1 else 0
          ) row
      ) arr;
    pref
  in
  let get_sum (x1, y1) (x2, y2) = 
    let _get (a, b) = try sum_pref.(a).(b) with _ -> 0 in 
    _get(x2, y2) - _get(x2, y1 - 1) - _get(x1 - 1, y2) + _get(x1 - 1, y1 - 1)
  in

  let check k =
    try
      Array.iteri ( fun x row -> 
        Array.iteri ( fun y el -> 
            let sum = try get_sum (x, y) (x + k, y + k) with _ -> 1 in
            if sum = 0 then raise True_exp
          ) row
      ) arr;
      false
    with
    | True_exp -> true
  in

  let p = ref (-1)
  and k = ref (1 + min n m) in
  while !k - !p > 1 do
    let midd = (!p + !k) / 2 in
    if check midd then
      p := midd
    else 
      k := midd
  done;

  !p + 1
;;


let bool_of_int x = not (x = 0)
;;
let bool_arr_of_int_arr = Array.map (bool_of_int)
;;
let bool_matrix_of_int_matrix = Array.map (bool_arr_of_int_arr)
;;

let tab = [|
  [|0; 1; 1; 1; 1; 1; 0; 0; 0; 1|];
  [|0; 1; 1; 1; 1; 1; 1; 0; 0; 1|];
  [|0; 1; 1; 1; 1; 1; 0; 0; 0; 1|];
  [|0; 1; 1; 1; 1; 1; 0; 0; 0; 0|];
  [|0; 1; 0; 0; 0; 1; 1; 0; 0; 0|];
  [|0; 0; 0; 1; 0; 0; 0; 1; 1; 1|];
  [|0; 0; 1; 1; 0; 0; 1; 0; 0; 1|];
  [|0; 1; 0; 1; 0; 1; 0; 1; 1; 1|];
  [|0; 0; 1; 1; 0; 1; 1; 1; 1; 1|];
  [|0; 0; 0; 1; 0; 1; 0; 1; 1; 1|];
|];;

let tab = bool_matrix_of_int_matrix tab
;;

let kw = kwadrat tab
;;

Printf.printf "Wynik: %d\n" kw
;;



