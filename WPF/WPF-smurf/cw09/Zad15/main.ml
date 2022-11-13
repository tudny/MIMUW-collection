
type kolor = Czerwona | Zielona

let kulki li =
  let cnt kol =
    List.fold_left
      (fun acc k -> acc + if k = kol then 1 else 0) 0
  in
  let cntZ = cnt Zielona li
  and cntC = cnt Czerwona li in
  let minCnt = min cntZ cntC in
  let maxK = if cntZ > cntC then Zielona else Czerwona in
  let jaki id =
    if id > minCnt * 2 then maxK
    else if id mod 2 = 1 then Zielona
    else Czerwona
  in
  let (zle, _) = List.fold_left
      (fun (acc, i) k ->
         if k <> jaki i then (acc + 1, i + 1)
         else (acc, i + 1))
    (0, 1)
    li
  in zle / 2
;;

kulki [Czerwona; Zielona; Zielona; Zielona; Zielona; Czerwona; Czerwona; Zielona];;


