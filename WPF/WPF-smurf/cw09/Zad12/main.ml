
let wybory lst =
  let trd (_, _, _, a) = a in
  let loop (ost, cnt, maks, maks_ele) el =
    let n_cnt = 1 + if ost = el then cnt else 0 in
    if n_cnt >= maks then
      (el, n_cnt, n_cnt, el)
    else
      (el, n_cnt, maks, maks_ele)
  in trd (List.fold_left (loop) (min_int, 1, 0, min_int) lst);;

let lista = [1; 2; 3; 3];;

wybory lista;;

let lista = [1; 1; 2; 3];;

wybory lista;;
