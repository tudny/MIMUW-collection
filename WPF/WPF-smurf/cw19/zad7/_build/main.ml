
type graph = int list array
;;

let polaczenie g u v w = 
  let n = Array.length g in 
  let inf = -1 in 
  let dis = Array.init 3 (fun _ -> Array.make n inf) in

  let calc_dis distance x = 
    let q = Queue.create () in
    Queue.add x q;
    distance.(x) <- 0;
    while not (Queue.is_empty q) do
      let akt = Queue.take q in
      List.iter (fun node ->
          if distance.(node) = inf then 
            distance.(node) <- distance.(akt) + 1;
            Queue.add node q
        ) g.(akt)
    done
  in

  let verts = [|u; v; w|] in
  Array.iter2 ( fun i dist -> calc_dis dist i) verts dis;

  let mini = ref max_int in
  Array.iteri (fun node _ -> 
      mini := min !mini (Array.fold_left ( fun acc dis -> acc + dis.(node) ) 0 dis)
    ) g; !mini
;;
