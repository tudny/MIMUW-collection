type 'a pri_queue =
    Node of 'a * 'a pri_queue * 'a pri_queue * int | Leaf

exception EmptyQueue

let empty_queue = Leaf

let is_empty q = q = Leaf

let size q =
  match q with
    Leaf -> 0 |
    Node (_, _, _, n) -> n

let get_root h =
  match h with
    Leaf -> raise EmptyQueue |
    Node (r, _, _, _) -> r

let set_root h r =
  match h with
    Leaf -> Node (r, Leaf, Leaf, 1) |
    Node (_, l, p, n) -> Node (r, l, p, n)

let rec put h x =
  match h with
    Leaf -> Node (x, Leaf, Leaf, 1) |
    Node (r, l, p, n) ->
    if(size l <= size p)
    then Node((max x r), (put l (min x r)), p, (n+1)) (* Mozna tak. *)
    else Node((max x r), l, (put p (min x r)), (n+1))

let getmax h =
  let rec del h =
    match h with
      Leaf -> h |
      (* Node (_, Leaf, Leaf, _) -> Leaf |  (* niepotrzebne *) *)
      Node (_, l, Leaf, _) -> l |
      Node (_, Leaf, p, _) -> p |
      Node (_, (Node (rl, _, _, _) as l), (Node (rp, _, _, _) as p), n) ->
      if rl >= rp then
        Node (rl, del l, p, n - 1)
      else
        Node (rp, l, del p, n - 1)
  in (get_root h, del h)



type graph = (int * int) list array
;;


let dijkstra (g : graph) (start : int) =
  let n = Array.length g in
  let distance = Array.make n max_int in
  let q = ref (put empty_queue (0, start)) in
  distance.(start) <- 0;
  while not (is_empty !q) do
    let ((dist, v), n_q) = getmax !q in
    q := n_q;
    if distance.(v) = -dist then begin
      List.iter ( fun (w, d) ->
          if distance.(w) > distance.(v) + d then begin
            distance.(w) <- distance.(v) + d;
            q := put !q (-distance.(w), w)
          end
        ) g.(v)
    end;
  done;
  distance
;;

