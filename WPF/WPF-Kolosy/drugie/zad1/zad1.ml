(* Na dole *)

type expr = Add of expr * expr | Var of string

let rec fold_expr f_add f_var = function
  | Add (l, r) -> f_add (fold_expr f_add f_var l) (fold_expr f_add f_var r)
  | Var z -> f_var z

let rec st = function Add(l, r) -> max (st l) (1 + st r) | Var _ -> 1


(*
let process (w_l, l) (w_p, p) =
  let maks_l_p = max w_l (w_p + 1)
  and maks_p_l = max w_p (w_l + 1) in
  if maks_l_p <= maks_p_l then
    (maks_l_p, Add (l, p))
  else
    (maks_p_l, Add (p, l))
;;

let process_end str =
  (1, Var(str))
;;
*)

(* Rozwiązanie | czasowo O(n) | pamięć O(n) | gdzie n - liczba wierzchołków drzewa *)
(* Dla każdego wierzchołka od dołu zastanawiam się *)
(* czy opłaca mi się zamienić synów miejscami *)
(* Dla dwóch możliwości liczę wynik *)
(* i w jego zależności zwracam odpowiednio stan normalny lub zamieniony *)
(* Minimalny stos otrzymamy, gry zawsze będziemy zwracać minimalny wybór *)
(* Drzewo budowane w pamięci zajmuje O(n) miejsca *)
(* Zmienne: (wl_l, l) - w_l wysokość stosu dla lewego drzewa, l - zbudowane lewe drzewo *)
(* (w_r, r) - analogicznie *)

let optymalizuj exp =
  let process (w_l, l) (w_p, p) =
    let maks_l_p = max w_l (w_p + 1)
    and maks_p_l = max w_p (w_l + 1) in
    if maks_l_p <= maks_p_l then
      (maks_l_p, Add (l, p))
    else
      (maks_p_l, Add (p, l))
  in
  let process_end str =
      (1, Var(str))
  in
  fold_expr (process) (process_end) exp |> snd
;;

(* ***** *)
(* TESTY *)
(* ***** *)

optymalizuj (Add (Var "a", Add (Var "b", Var "c")));;

let ex = (Add (Var "a", Add(Var "b", Add(Var "c", Add(Var "d", Var "e")))));;
let ex_opt = optymalizuj ex;;
st ex;;
st ex_opt;;

