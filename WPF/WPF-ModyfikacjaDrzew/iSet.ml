(* Autor: Aleksander Tudruj *)
(* Review: choroba *)

(* typ odpowiadający za przedział całkowitoliczbowy *)
type interval = int * int

(* Empty *)
(* Node (left_tree, (a, b), right_tree, height, number_of_elements_int_tree) *)
type t =
  | Empty
  | Node of t * interval * t * int * int;;

(* funkcja zwracająca wysokość drzewa *)
let height = function
  | Empty -> 0
  | Node(_, _, _, h, _) -> h

(* funkcja zwracająca wielkość drzewa *)
let size = function
  | Empty -> 0
  | Node(_, _, _, _, n) -> n

(* własny operator dodawania uwzględniający overflow *)
let (+!) a b =
  if a + b < 0 then max_int else a + b

(* operator porównujący "ostrą nierówność", która zachodzi gdy a + 1 < b *)
let (<<) a b =
  if a > 0 then a +! 1 < b
  else a + 1 < b

(* długość przedziału a b *)
let length (a, b) =
  if b - a + 1 <= 0 then max_int else b - a + 1

(* sprawdzenie czy liczba jest w przedziale *)
let in_range x (a, b) =
  a <= x && x <= b

(* proste tworzenie drzewa z wierzchołka i dwóch poddrzew *)
let make l i r =
  Node (l, i, r, max (height l) (height r) + 1,
        (size l) +! (size r) +! (length i))

(* sprawdzenie czy drzewo jest puste *)
let is_empty s =
  s = Empty

(* zwraca puste drzewo *)
let empty = Empty

(* funkcja balansująca (jak w oryginalnej bibliotece) *)
(* O(1) *)
let bal l k r =
  let hl = height l in
  let hr = height r in
  if hl > hr + 2 then
    match l with
    | Node (ll, lk, lr, _, _) ->
      if height ll >= height lr then make ll lk (make lr k r)
      else
        (match lr with
         | Node (lrl, lrk, lrr, _, _) ->
           make (make ll lk lrl) lrk (make lrr k r)
         | Empty -> assert false)
    | Empty -> assert false
  else if hr > hl + 2 then
    match r with
    | Node (rl, rk, rr, _, _) ->
      if height rr >= height rl then make (make l k rl) rk rr
      else
        (match rl with
         | Node (rll, rlk, rlr, _, _) ->
           make (make l k rll) rlk (make rlr rk rr)
         | Empty -> assert false)
    | Empty -> assert false
  else make l k r

(* minimalny element drzewa *)
(* O(log n), n - wielkość drzewa *)
let rec min_elt = function
  | Node (Empty, k, _, _, _) -> k
  | Node (l, _, _, _, _) -> min_elt l
  | Empty -> raise Not_found

(* maksymalny element drzewa *)
(* O(log n), n - wielkość drzewa *)
let rec max_elt = function
  | Node (_, k, Empty, _, _) -> k
  | Node (_, _, r, _, _) -> max_elt r
  | Empty -> raise Not_found

(* usuwa minimalny element drzewa *)
(* O(log n), n - wielkość drzewa *)
let rec remove_min_elt = function
  | Node (Empty, _, r, _, _) -> r
  | Node (l, k, r, _, _) -> bal (remove_min_elt l) k r
  | Empty -> invalid_arg "ISet.remove_min_elt"

(* usuwa maksymalny element drzewa *)
(* O(log n), n - wielkość drzewa *)
let rec remove_max_elt = function
  | Node (l, _, Empty, _, _) -> l
  | Node (l, k, r, _, _) -> bal l k (remove_max_elt r)
  | Empty -> invalid_arg "ISet.remove_max_elt"

(* proste dodanie elementu zakładające, że nie możne go połączyć *)
(* z żadnym istniejącym w tym drzewie elemencie *)
(* O(log n) *)
let rec add_one_simple (a, b) = function
  | Empty -> make Empty (a, b) Empty
  | Node (l, (x, y), r, h, _) ->
    if b << x then
      bal (add_one_simple (a, b) l) (x, y) r
    else if y << a then
      bal l (x, y) (add_one_simple (a, b) r)
    else
      invalid_arg "ISet.add_one_simple"

(* łączenie dwóch niezbalansowanych drzew i wierzchołka *)
(* amortyzuje się do O(log n) *)
let rec join l i r =
  match (l, r) with
  | (Empty, _) -> add_one_simple i r
  | (_, Empty) -> add_one_simple i l
  | (Node(ll, lv, lr, lh, _), Node(rl, rv, rr, rh, _)) ->
    if lh > rh + 2 then bal ll lv (join lr i r) else
    if rh > lh + 2 then bal (join l i rl) rv rr else
      make l i r

(* funkcja rozdzielająca wedle polecenia *)
(* O(log n) *)
let split a s =
  let rec loop = function
    | Empty -> (Empty, false, Empty)
    | Node (l, (x, y), r, _, _) ->
      if a < x then
        let (l_subtree, is, r_subtree) = loop l in
        (l_subtree, is, join r_subtree (x, y) r)
      else if y < a then
        let (l_subtree, is, r_subtree) = loop r in
        (join l (x, y) l_subtree, is, r_subtree)
      else
        let l_st =
          if x < a then add_one_simple (x, a - 1) l else l
        and r_st =
          if a < y then add_one_simple (a + 1, y) r else r
        in (l_st, true, r_st)
  in loop s

(* dodawanie elementu do drzewa *)
(* amortyzuje się do O(log n) *)
let add (a, b) t =
  let (l_st, _, _) = split a t
  and (_, _, r_st) = split b t in
  match l_st, r_st with
  | Empty, Empty -> make Empty (a, b) Empty
  | (Empty, Node(_)) ->
    let (x, y) = min_elt r_st in
    if b << x then
      join empty (a, b) r_st
    else
      join empty (a, y) (remove_min_elt r_st)
  | (Node(_), Empty) ->
    let (x, y) = max_elt l_st in
    if y << a then
      join l_st (a, b) empty
    else
      join (remove_max_elt l_st) (x, b) empty
  | (Node(_), Node(_)) ->
    let (ll, min_i) =
      let (x, y) = max_elt l_st in
      if y << a then (l_st, a) else (remove_max_elt l_st, x)
    and (rr, max_i) =
      let (x, y) = min_elt r_st in
      if b << x then (r_st, b) else (remove_min_elt r_st, y)
    in join ll (min_i, max_i) rr

(* usuwanie elementu z drzewa *)
(* O(log n) *)
let remove (a, b) s =
  let (l_st, _, _) = split a s
  and (_, _, r_st) = split b s in
  match l_st, r_st with
  | Empty, Empty -> Empty
  | Empty, Node (_) ->
    let i = min_elt r_st in
    join l_st i (remove_min_elt r_st)
  | Node (_), _ ->
    let i = max_elt l_st in
    join (remove_max_elt l_st) i r_st

(* sprawdzenie czy element jest w drzewie *)
(* O(log n) *)
let rec mem x s =
  let rec loop = function
    | Empty -> false
    | Node (l, (a, b), r, _, _) ->
      in_range x (a, b) ||
      if x < a then loop l else loop r
  in loop s

(* funkcja z polecenia *)
(* O(n) *)
let rec iter f s =
  let rec loop = function
    | Empty -> ()
    | Node (l, i, r, _, _) ->
      loop l; f i; loop r
  in loop s

(* funkcja z polecenia *)
(* O(n) *)
let rec fold f s a =
  let rec loop acc = function
    | Empty -> acc
    | Node (l, i, r, _, _) ->
      loop (f i (loop acc l)) r in
  loop a s

(* funkcja z polecenia *)
(* O(n) *)
let elements s =
  let rec loop acc = function
    | Empty -> acc
    | Node (l, i, r, _, _) ->
      loop (i :: loop acc r) l in
  loop [] s

(* funkcja z polecenia *)
(* O(log n) *)
let below x s =
  let rec loop acc = function
    | Empty -> acc
    | Node (l, (a, b), r, _, _) ->
      if x < a then loop acc l
      else if in_range x (a, b) then (size l) +! length (a, x) +! acc
      else loop ((size l) +! length (a, b) +! acc) r
  in loop 0 s
