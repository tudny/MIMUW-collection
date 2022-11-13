
(* Czas: O(n) *)
(* Pamięc: O(n) *)

(* Brutalne zapisanie na tablicy, odwracanie fragmentów i przemięcie wskaźników *)

(* Zadanie 1 *)

type 'a option = None | Some of 'a
type 'a elem = {v: 'a; mutable prev: 'a lista; mutable next: 'a lista}
and 'a lista = 'a elem option;;

let get_some = function
  | Some x -> x
  | None -> invalid_arg "noneeee"
;;

let poodwracaj li leny =
  let n = List.fold_left (+) 0 leny in
  let temp = Array.make n None in

  let x = ref li in
  let cnt = ref 0 in
  while !x <> None do
    temp.(!cnt) <- !x;
    x := (get_some !x).next;
    cnt := !cnt + 1
  done;

  let rev x y arr =
    let x = ref x
    and y = ref y in
    while !x < !y do
      let temp = arr.(!x) in
      arr.(!x) <- arr.(!y);
      arr.(!y) <- temp;
      x := !x + 1;
      y := !y - 1;
    done
  in

  let w = ref 0 in
  List.iter ( fun len ->
      let v = !w + len - 1 in
      rev !w v temp;
      w := !w + len
    ) leny;

  Array.iteri ( fun i el ->
      (get_some el).next <- if i = n - 1 then None else temp.(i + 1);
      (get_some el).prev <- if i = 0 then None else temp.(i - 1)
    ) temp;

  temp.(0)
;;


(* testy *)


let rec d1 = Some {v = 4; prev = None; next = d2}
  and d2 = Some {v = 5; prev = d1; next = d3}
  and d3 = Some {v = 6; prev = d2; next = d4}
  and d4 = Some {v = 7; prev = d3; next = d5}
  and d5 = Some {v = 8; prev = d4; next = d6}
  and d6 = Some {v = 9; prev = d5; next = None};;
let d = d1;;

poodwracaj d [3; 2; 1];;

d;;




