
(* czas: O(n * m * log^{*}(n*m)) *)
(* pamięć: O(n * m) *)

(* log^{*} to złożoność Find and Union *)

(* Idę od tyłu po kolei po ścianach *)
(* Na poczatku zakładam, że wszystkie są postawione *)
(* i jakby burzę ściany *)
(* Jeżeli zburzenie ściany sprawiło, że dwa pola się nie połączyły *)
(* To była to interesująca nas sciana *)

(* Zadanie 2 *)

let labirynt n m sciany =
  let rep = Array.make_matrix n m (0, 0) in
  let size = Array.make_matrix n m 1 in
  let rec find (x, y) =
    if (x, y) = rep.(x).(y) then (x, y) else begin
      rep.(x).(y) <- find rep.(x).(y);
      rep.(x).(y)
    end
  in
  for x = 0 to n - 1 do
    for y = 0 to m - 1 do
      rep.(x).(y) <- (x, y);
      size.(x).(y) <- 1
    done
  done;
  let union (x, y) (a, b) =
    let (a, b) = find (a, b)
    and (x, y) = find (x, y) in
    if (a, b) = (x, y) then false else begin
      if size.(a).(b) > size.(x).(y) then begin
        rep.(x).(y) <- (a, b);
        size.(a).(b) <- size.(a).(b) + size.(x).(y)
      end else begin
        rep.(a).(b) <- (x, y);
        size.(x).(y) <- size.(x).(y) + size.(a).(b)
      end;
      true
    end
  in
  let sciany = List.rev sciany in
  let res = ref [] in
  List.iter ( fun (x, y, v) ->
      let ((a, b), (c, d)) =
        if v then
          (find (x, y), find (x, y - 1))
        else
          (find (x, y), find (x - 1, y)) in

      if not (union (a, b) (c, d)) then
        res := (x, y, v) :: !res
    ) sciany;
  !res |> List.rev
;;


(* Test *)

let test = labirynt 3 2 [(2,1,false); (2,1,true); (2,0,false); (1,1,true);
                         (0,1,true); (1,1,false); (1,0,false) ] in

assert (test = [(1,1,true); (2,1,false)])

























