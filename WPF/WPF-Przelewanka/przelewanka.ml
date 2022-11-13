(* ************************** *)
(*        Przelewanka         *)
(* Autor: Aleksander Tudruj   *)
(* Review: Mikołaj Piróg      *)

(* Wyjątek podnoszony w razie znalezienia wyniku. *)
(* Skraca czas wykonania programu. *)
exception Found of int

(* Rozpatruję wszystkie możliwoście przelewania/wlewania/wylewania  *)
(* i podobnie do działania algorytmu BFS szukam najkrótszej ścieżki *)
(* w grafie. *)
let przelewanka cups =
  let cups_nr = Array.length cups
  in
  if cups_nr = 0 then 0 else

  let array_nwd =
    let rec nwd x y =
      if y = 0 then x
      else nwd y (x mod y)
    in Array.fold_left nwd 0
  in

  let array_split arr =
    let len = Array.length arr in
    let part which = Array.init len (fun i -> which arr.(i)) in
    (part fst, part snd)
  in

  let heights, result = array_split cups
  in

  let d = array_nwd heights in
  if d = 0 then 0 else

  (* Sprawdzam warunki na istnienie wyniku. *)
  let can_be_done =
    Array.for_all (fun h -> h mod d = 0) result &&
    Array.exists (fun (h, e) -> h = e || e = 0) cups
  in

  if not (can_be_done) then -1 else

  (* Symulacja polegająca na dolaniu wody do pełna do kubka [x] *)
  let until_full x state =
    let state = Array.copy state in
    state.(x) <- heights.(x);
    state
  in

  (* Symulacja polegająca na wylaniu wody z kubka [x] *)
  let until_empty x state =
    let state = Array.copy state in
    state.(x) <- 0;
    state
  in

  (* Symulacja polegająca przelaniu wody z kubka [x] do kubka [y] *)
  let pour state_def y x =
    let state = Array.copy state_def in
    let left = heights.(y) - state.(y) in
    if left >= state.(x) then begin
      state.(y) <- state.(y) + state.(x);
      state.(x) <- 0;
    end else begin
      state.(y) <- heights.(y);
      state.(x) <- state.(x) - left
    end;
    state
  in

  (* Kolejka niczym z algorytmu BFS; tak samo mapa odwiedzonych wierzchołków *)
  let q = Queue.create () in
  (* O dziwo to magiczne ustawienie Hashtbl daje najlepsze czasy *)
  let vis = Hashtbl.create (cups_nr * 202_020) ~random:true in
  let start_state = Array.make cups_nr 0 in

  (* Rozpatrując nową drogę w przelewaniu wody, sprawdzam czy nie osiągnąłem *)
  (* jej już wczesniej. Jeżeli nie, to ustawiam na odwiedzone (a konkretnie  *)
  (* ustawiam dystans do tego wierzchołka od stanu poczatkowego) i dodaję    *)
  (* do kolejki.                                                             *)
  let process new_way dis =
    if new_way = result then raise (Found dis);
    if not (Hashtbl.mem vis new_way) then begin
      Hashtbl.add vis new_way dis;
      Queue.add new_way q
    end
  in

  try
    process start_state (0);

    while not (Queue.is_empty q) do
      let front = Queue.take q in
      let dis = Hashtbl.find vis front in

      (* Sprawdzam wszystkie możliwości.                       *)
      (* Nie rozpatruję tych, które nie mają logicznego sensu. *)
      Array.iteri ( fun ele height ->
          (* Wylanie *)
          if height <> 0 then
            process (until_empty ele front) (dis + 1);

          (* Dolanie *)
          if height <> heights.(ele) then
            process (until_full ele front) (dis + 1);

          (* Przelanie *)
          Array.iteri ( fun ele2 _ ->
              if ele <> ele2 && heights.(ele2) <> 0 &&
               heights.(ele) <> height then begin
                process (pour front ele ele2) (dis + 1)
              end
            ) front
        ) front
    done;
    -1
  with
  | Found x -> x

