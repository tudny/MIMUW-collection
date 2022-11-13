
let diagonal arr =

  let n = Array.length arr in

  let calc_diag (a, b) = sqrt (a *. a +. b *. b) in

  let (f, s) = arr.(0) in
  if f -. s <= 0. then calc_diag arr.(0) else

  let (f, s) = arr.(n - 1) in
  if f -. s >= 0. then calc_diag arr.(n - 1) else

  let p = ref 0
  and k = ref (n - 1) in

  while !p + 1 < !k do
    let midd = (!p + !k) / 2 in
    let (f, s) = arr.(midd) in
    if f -. s >= 0. then
      p := midd
    else
      k := midd
  done;
  let diag1 = calc_diag arr.(!p) in

  let p = ref 0
  and k = ref (n - 1) in

  while !p + 1 < !k do
    let midd = (!p + !k) / 2 in
    let (f, s) = arr.(midd) in
    if f -. s > 0. then
      p := midd
    else
      k := midd
  done;
  let diag2 = calc_diag arr.(!k) in
  max diag1 diag2
;;


let diag =
  diagonal [|(24., 1.); (12., 2.); (6., 4.); (2., 12.); (1., 24.)|];;

Printf.printf "Diag: %f\n" diag;;
