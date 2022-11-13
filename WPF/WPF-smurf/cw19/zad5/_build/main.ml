
exception Brzeg
;;

let sadzawka mapa jas = 
  let (a, b) = jas in
  let n = Array.length mapa
  and m = Array.length mapa.(0) in

  let in_range x (a, b) = 
    a <= x && x < b
  in

  let check moment = 

    let vis = Array.make_matrix n m false in

    let rec _dfs (x, y) = 
      vis.(x).(y) <- false;
      let dx = [0; 1; 0; -1]
      and dy = [1; 0; -1; 0] in
      List.iter2 ( fun dx dy -> 
          let nx = x + dx
          and ny = y + dy in 
          if in_range nx (0, n) && in_range ny (0, m) then begin
            if not vis.(nx).(ny) && mapa.(nx).(ny) >= moment then
              _dfs (nx, ny)
          end else
            raise Brzeg 
        ) dx dy
    in

    try 
      _dfs (a, b); false
    with
    | Brzeg -> true
  in
  
  let p = ref 0
  and k = ref (mapa.(a).(b)) in

  while !k - !p > 1 do
    let midd = (!k + !p) / 2 in
    if check midd then
      p := midd
    else 
      k := midd
    done;
    !p
;;


let lod = [|
  [|4; 5;  1; 2|];
  [|7; 4; 10; 3|];
  [|1; 3;  2; 1|];
  [|2; 2;  2; 2|]
|];;

(* let lod = [|
  [|2; 2; 2|];
  [|2; 10; 2|];
  [|2; 2; 2|];
|];; *)

let jas = (1, 1);;

let res = sadzawka lod jas;;

Printf.printf "sadzawka: %d\n" res
;;
