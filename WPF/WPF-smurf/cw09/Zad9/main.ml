
let zlecenia zl =
  let loop (moment_skonczenia_ostatnego, acc) (pocz, ile) =
    let akt_pocz = if pocz > moment_skonczenia_ostatnego then
        pocz else moment_skonczenia_ostatnego in
    let akt_kon = akt_pocz + ile in
    (akt_kon, akt_kon - pocz :: acc) in
  List.rev (snd (List.fold_left (loop) (min_int, []) zl));;

zlecenia [(-1, 1); (2, 2); (3, 3); (4, 2); (10, 2)];;
assert (zlecenia [(-1, 1); (2, 2); (3, 3); (4, 2); (10, 2)] = [1; 2; 4; 5; 2]);;
