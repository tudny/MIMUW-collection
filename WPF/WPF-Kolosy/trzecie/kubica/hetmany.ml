(* Problem 8 hetmanów *)
#use "backtracking.ml";;
open List;;

module Hetman = 
  struct
    (* Pozycje hetmanów na szachownicy. *)
    type result = (int * int) list

    (* Typ konfiguracji na szachownicy.                      *)
    (* Trójka: wielko¶æ planszy, liczba pocz±tkowych kolumn, *)
    (* w których nale¿y postawiæ hetmany, lista ich pozycji. *)
    type config = int * int * result

    (* Pusta plansza rozmiaru n. *)
    let empty n = (n, n, [])

    (* Pozycje hetmanów w konfiguracji. *)
    let extract (_, _, l) = l

    (* Czy gotowe rozwi±zanie. *)
    let final (_, k, _) = k=0

    (* Funkcja okre¶laj±ca, czy dwa hetmany siê szachuj± *)
    let szach (x1, y1) (x2, y2) = 
      x1=x2 or y1=y2 or x1+y1=x2+y2 or x1-y1=x2-y2

    (* Czy mo¿na dostawiæ hetmana h do ustawionych ju¿ hetmanów het *)
    let mozna_dostawic h het = 
      fold_left (fun ok x -> ok && not (szach h x)) true het
	  
    (* Lista kolejnych liczb ca³kowitych od-do. *)
    let rec ints a b = if a > b then [] else a :: ints (a+1) b;;

    (* Konfiguracje powsta³e przez dostawienie hetmana w kolumnie k. *)
    let iter p (n, k, l) = 
      let r = filter (fun i -> mozna_dostawic (k, i) l) (ints 1 n)
      in List.iter (fun i -> p (n, k-1, (k,i)::l)) r
  end;;

(* Instancja backtrackingu znajduj±ca jedno rozwi±zanie. *)
module H = Bt_Solver (Hetman);;

(* Procedura znajduj±ca ustawienie hetmanów. *)
let hetman n = H.onesolution (Hetman.empty n);;

(* Procedura znajduj±ca ustawienie hetmanów. *)
let hetmany n = H.solutions (Hetman.empty n);;

(* Testy. *)
assert (hetman 0 = []);;
assert (hetman 1 = [(1, 1)]);;
assert (let w = hetman 4 
        in  w = [(1, 3); (2, 1); (3, 4); (4, 2)] || 
            w = [(1, 2); (2, 4); (3, 1); (4, 3)]);;

assert (hetmany 0 = [[]]);;
assert (hetmany 1 = [[(1, 1)]]);;
assert (let w = hetmany 4 
        in  w = [[(1, 3); (2, 1); (3, 4); (4, 2)]; [(1, 2); (2, 4); (3, 1); (4, 3)]] || 
            w = [[(1, 2); (2, 4); (3, 1); (4, 3)]; [(1, 3); (2, 1); (3, 4); (4, 2)]]);;



