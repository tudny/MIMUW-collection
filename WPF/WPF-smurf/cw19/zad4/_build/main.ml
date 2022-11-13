
(* Bin Search po wyniku *)
let wysokosc mapa = 
  (* Mapa nie może być pusta *)
  let n = Array.length mapa in
  let m = Array.length mapa.(0) in

  let remove_duplicates li = 
    match li with
    | [] | [_] -> li 
    | h :: t -> 
      List.fold_left (fun (last, acc) ele -> 
        if ele <> last then
          (ele, ele :: acc)
        else
          (ele, acc)
        ) (h, [h]) t |> snd |> List.rev
  in

  let all_heights =
    List.flatten (Array.to_list (Array.map (Array.to_list) mapa)) |>
    List.sort compare |>
    remove_duplicates |> 
    Array.of_list in
  let h = Array.length all_heights in

  let check maximal = 
    if mapa.(0).(0) > maximal then false else
    let vis = Array.make_matrix n m false in
    let rec dfs_max (x, y) = 
      vis.(x).(y) <- true;
      let dx_t = [0; 1; 0; -1]
      and dy_t = [1; 0; -1; 0] in
      List.iter2 (fun dx dy -> 
          let nx = x + dx
          and ny = y + dy in 
          if 0 <= nx && nx < n && 0 <= ny && ny < m && 
                not vis.(nx).(ny) && mapa.(nx).(ny) <= maximal then
            dfs_max (nx, ny)
        ) dx_t dy_t in 
    dfs_max (0, 0);
    vis.(n - 1).(m - 1)
  in

  let p = ref 0
  and k = ref (h - 1) in

  while !k - !p > 1 do
    let mid = (!k + !p) / 2 in 
    let max_h = all_heights.(mid) in
    if check max_h then
      k := mid
    else
      p := mid 
  done;
  all_heights.(!k)
;;

let m = 
[|
  [|2; 1; 7; 3|];
  [|3; 1; 4; 1|];
  [|8; 1; 3; 2|];
  [|2; 7; 4; 1|]
|];;

Printf.printf "wys: %d\n" (wysokosc m);;
(* wys: 3 *)

let m = 
  [|
    [|2; 7; 4; 1|];
    [|8; 1; 3; 2|];
    [|3; 1; 4; 1|];
    [|2; 1; 7; 3|]
  |];;

Printf.printf "wys: %d\n" (wysokosc m);;
(* wys: 7 *)
