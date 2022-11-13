
type graf = int list array
;;

let dfs (graph : graf) (start : int) =
  let n = Array.length graph in
  let visited = Array.make n false in

  let rec _dfs v =
    visited.(v) <- true;
    List.iter (fun w ->
        if not visited.(w) then
          _dfs w
      ) graph.(v)
  in
  _dfs start
;;


let components (graph : graf) =
  let n = Array.length graph in
  let visited = Array.make n false in
  let num_comp = ref 0 in

  let rec _dfs v =
    visited.(v) <- true;
    List.iter (fun w ->
        if not visited.(w) then
          _dfs w
      ) graph.(v)
  in

  for i = 0 to (n - 1) do
    if not visited.(i) then begin
      num_comp := !num_comp + 1;
      _dfs i
    end
  done;

  !num_comp;;


let g =
  [|
    [1; 2]; (* 0 *)
    [0];    (* 1 *)
    [0];    (* 2 *)
    [4];    (* 3 *)
    [3];    (* 4 *)
    [5; 6]; (* 5 *)
    [5]     (* 6 *)
  |];;
(* 2 - 0 - 1    3 - 4    5 - 6 *)

components g;;

