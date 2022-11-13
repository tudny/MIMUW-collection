
type graph = int list array
;;

let bfs (g : graph) (start : int) =
  let n = Array.length g in
  let visited = Array.make n false in
  let q = Queue.create () in
  visited.(start) <- true;
  Queue.push start q;
  while not (Queue.is_empty q) do
    let v = Queue.take q in
    List.iter (fun w ->
        if not visited.(w) then begin
          Queue.push w q;
          visited.(w) <- true;
        end
      ) g.(v)
  done
;;

let bfs_odleglosci (g : graph) (start : int) =
  let n = Array.length g in
  let distance = Array.make n (-1) in
  let q = Queue.create () in
  distance.(start) <- 0;
  Queue.push start q;
  while not (Queue.is_empty q) do
    let v = Queue.take q in
    List.iter (fun w ->
        if distance.(w) = (-1) then begin
          Queue.push w q;
          distance.(w) <- distance.(w) + 1;
        end
      ) g.(v)
  done;
  distance
;;


