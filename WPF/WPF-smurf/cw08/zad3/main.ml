
let odwrotnosc f x =
  let eps = 1e-6 in
  let rec bin_s b e =
    let m = (b +. e) /. 2. in
    let value = f m in
    if value -. x < eps then m
    else if value < x then bin_s m e
    else bin_s b m
  in if x > 0. then bin_s 0. x else bin_s x 0.
