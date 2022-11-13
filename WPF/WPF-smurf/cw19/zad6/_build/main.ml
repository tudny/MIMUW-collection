
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




let wyspa (mapa : bool array array) =
  let n = Array.length mapa in
  let m = Array.length mapa.(0) in
  let neg f x = not (f x) in
  let int_of_bool x = if x then 1 else 0 in
  let poza (x, y) = x = 0 || x = n + 1 || y = 0 || y = m + 1 in
  let czy_woda (x, y) =
    poza (x, y) || not mapa.(x).(y) in
  let czy_lad = neg czy_woda in
  let in_range (x, y) = x >= 0 && x <= n + 1 && y >= 0 && y <= m + 1 in
  let warstwa = Array.make_matrix (n + 2) (m + 2) (-1) in
  let set_war (x, y) nr = warstwa.(x + 1).(y + 1) <- nr in
  let get_war (x, y) = warstwa.(x + 1).(y + 1) in
  set_war (0, 0) 0;
  (* (nr_warstwy, lad0/woda1, (x, y)) *)
  let q = ref (put empty_queue (0, 1, (0, 0))) in
  while not (is_empty !q) do
    let ((nr_war, typ, (x, y)), n_q) = getmax !q in
    Printf.printf "(%d, %d)\n" x y;
    q := n_q;
    let dx = [0; 1; 0; -1]
    and dy = [1; 0; -1; 0] in
    List.iter2 ( fun dx dy ->
        let np = (x + dx, y + dy) in
        Printf.printf "n(%d, %d)\n" (fst np) (snd np);
        if in_range np then
        Printf.printf "n2(%d, %d)\n" (fst np) (snd np);
        if get_war np = (-1) then begin
          if czy_lad (x, y) && czy_lad np then
            set_war np (get_war (x, y))
          else if czy_lad (x, y) && czy_woda np then
            set_war np (get_war (x, y) + 1)
          else if czy_woda (x, y) && czy_woda np then
            set_war np (get_war (x, y))
          else if czy_woda (x, y) && czy_lad np then
            set_war np (get_war (x, y));

          q := put !q (get_war np, int_of_bool (czy_woda np), np)
        end
      ) dx dy
  done;
  Array.fold_left (Array.fold_left (max)) 0 warstwa
;;


let wys = [|
  [|0; 0; 0; 0; 0; 0|];
  [|0; 1; 1; 1; 1; 0|];
  [|1; 1; 0; 0; 1; 0|];
  [|1; 0; 1; 0; 1; 0|];
  [|0; 1; 0; 1; 1; 0|];
  [|0; 0; 1; 1; 0; 0|];
  [|0; 0; 0; 0; 0; 0|];
|]

let bool_of_int x = x <> 0;;
let wys = Array.map (Array.map (bool_of_int)) wys;;

Printf.printf "%d\n" (wyspa wys);;

