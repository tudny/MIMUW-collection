
(* Graf na listach sąsiedztwa (node, distance) *)
type graph_lists = (int * int) list array;;

(* Graf na tablicy sąsiedztwa *)
type graph = int array array;;
(* N x N <- rozmiar *)
(* Tablica krawędzi *)
(* a.(x).(y) = *)
(* <int> jeżeli istnieje krawędź x->y lub x=y *)
(* max_int wpp *)

let convert_to_graph (graf : graph_lists) : graph =
  let n = Array.length graf in
  let g = Array.make_matrix n n max_int in
  
  Array.iteri ( fun i lst -> 
      List.iter ( fun (j, d) ->
        g.(i).(j) <- min g.(i).(j) d 
      ) lst;
      g.(i).(i) <- 0
    ) graf;
    g
;;

let floyd_warshall (graf : graph) : (int array array) =
  let n = Array.length graf in  
  let dis = Array.make_matrix n n (1_000_000_000) in 
  Array.iteri (fun i row -> 
      Array.iteri (fun j el -> 
          if el <> max_int then
            dis.(i).(j) <- graf.(i).(j)
        ) row
    ) graf;
  
  for k = 0 to n - 1 do
    for i = 0 to n - 1 do
      for j = 0 to n - 1 do
        if dis.(i).(j) > dis.(i).(k) + dis.(k).(j) then
          dis.(i).(j) <- dis.(i).(k) + dis.(k).(j)
      done
    done
  done;
  dis
;;


