(*
  Dana jest prostokątna mapa wysokości górzystego terenu, w postaci prostokątnej
  tablicy dodatnich liczb całkowitych, o wymiarach N×M.
  Chcemy przejść z pola o współrzędnych(0,0) na pole o współrzędnych(N−1,M−1),
  ale nie chcemy wspinać się zbyt wysoko.
  (Możemy się przesuwać w kierunkach N, W, S, E.)
  Napisz procedurę wysokosc:int array array→int,
  która dla danej mapy terenu określi
  minimalną największą wysokość, na którą musimy wejść w trakcie podróży.
*)

type graf = int array array;;

let wysokosc (mapa : graf) =
  let n = Array.length mapa in
  let m = Array.length mapa.(0) in

  let wysokosci = Array.make ((n * m) + 2) 0 in

  (* inicjacja pomocniczej listy z wysokosciami *)
  let pom_i = ref 1 in
  for i = 0 to n - 1 do
    for j = 0 to m - 1 do
      wysokosci.(!pom_i) <- mapa.(i).(j);
      pom_i := !pom_i + 1
    done;
  done;
  Array.sort compare wysokosci;

  let warunek k =
    let visited = Array.make_matrix n m false in
    let rec _dfs (x, y) =
      if mapa.(y).(x) > k then () else
        begin
          visited.(x).(y) <- true;
          if (x > 0) && (visited.(y).(x - 1)) then _dfs (x - 1, y);
          if (x < m - 1) && (visited.(y).(x + 1)) then _dfs (x + 1, y);
          if (y > 0) && (visited.(y - 1).(x)) then _dfs (x, y - 1);
          if (y < n - 1) && (visited.(y + 1).(x)) then _dfs (x, y + 1);
        end
    in _dfs (0, 0); visited.(n - 1).(m - 1)
  in

  let p = ref 0 and k = ref ((n * m) + 1) in

  while !k - !p > 1 do
    let mid = wysokosci.((!k + !p) / 2) in

    if warunek mid then
      k := mid
    else
      p := mid;
  done;
  !k;;

wysokosc
[|
    [| 1; 2; 3 |];
    [| 1; 2; 3 |];
    [| 1; 1; 1 |];
    [| 1; 3; 1 |];
    [| 2; 1; 1 |]
|];;