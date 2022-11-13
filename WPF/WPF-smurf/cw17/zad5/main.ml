
let skarbonki arr =
  let n = Array.length arr
  in
  let vis = Array.make n false
  in
  let graph = Array.init n (fun _ -> [])
  in
  Array.iteri ( fun i ele ->
      graph.(i) <- ele :: graph.(i);
      graph.(ele) <- i :: graph.(i)
    ) arr;

  let rec dfs x =
    vis.(x) <- true;
    List.iter (fun ele ->
        if not vis.(ele) then
          dfs ele
      ) graph.(x)
  in

  let cnt = ref 0
  in
  let incr x =
    x := !x + 1
  in
  Array.iteri (fun i _ ->
      if not vis.(i) then begin
        incr cnt;
        dfs i
      end
    ) graph;
  !cnt
;;



skarbonki [|1; 0; 1; 3|];;

