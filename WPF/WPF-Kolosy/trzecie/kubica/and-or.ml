(* Drzewa AND-OR. *)

(* Przydatne operacje na listach. *)
open List;;
let infty = max_int;;
let sum l = fold_left (+) 0 l;;
let minl l = fold_left min (hd l) (tl l);;

(* Mapy. *)
#use "mapy.ml";;

(* Gracze *)
type player = And | Or;;
let i = Or;;

(* Przeciwnik gracza who *)
let other = function And -> Or | Or -> And;;

(* Trójwarto¶ciowa logika: prawda, fa³sz i warto¶æ nieokre¶lona. *)
type boolean = True | False | Undefined;;

let and3 x y = 
  match (x, y) with 
    (False, _) -> False |
    (_, False) -> False |
    (True, _)  -> y |
    (_, True)  -> x |
    _          -> Undefined;;

let or3 x y = 
  match (x, y) with 
    (True, _)  -> True |
    (_, True)  -> True |
    (False, _) -> y |
    (_, False) -> x |
    _          -> Undefined;;

(* Wyniki gier. *)
let won = True;;
let lost = False;;

(* Gra w kamienie - czê¶æ zale¿na od gry *)
#use "gra-w-kamienie.ml";;

    
(* Drzewo AND-OR jest reprezentowane za pomoc± mapy przyporz±dkowuj±cej *)
(* parom (sytuacja, gracz) n-tki postaci:                               *)
(* - warto¶æ: True / False / Undefined                                  *)
(* - minimalna liczba el. falsyfikuj±cych,                              *)
(* - minimalna liczba el. potwierdzaj±cych,                             *)
(* - dla rozwiniêtych: lista trójek (ruch, sytuacja, gracz).            *)


(* Selektory *)
let vertex dag s p = apply dag (s, p);;

let value (v, _, _, _) = v;;
let min_f (_, f, _, _) = f;;
let min_t (_, _, t, _) = t;;
let succ  (_, _, _, s) = s;;

let succ_move  (m, _, _)   = m;;
let succ_board (_, b, _)   = b;;
let succ_who   (_, _, who) = who;;
let follow dag x = vertex dag (succ_board x) (succ_who x);;

(* Wêze³ jest rozwiniêty, je¿eli ma synów lub ma warto¶æ. *)
let is_expanded v = (value v <> Undefined) || (succ v <> []);;


(* Operacje na DAG-u *)

(* Czy ju¿ jest taki wierzcho³ek? *)
let is_present dag s who = dom dag (s, who);;

(* Dodanie do DAG-u wêz³a o parametrach:         *)
(* - s   - sytuacja,                             *)
(* - who - gracz,                                *)
(* - v   - wartosc,                              *)
(* - e   - lista trójek (ruch, sytuacja, gracz). *)
let insert dag s who v e =

  (* Rozmiar najmniejszego zbioru falsyfikuj±cego/potwierdzaj±cego *)
  (* (w zale¿no¶ci od parametrów selector/cumulator) dla synów.    *)
  let minft selector cumulator =
    (*  Lista liczb do przeliczenia, pochodz±cych od  *)
    (* synów z nieokre¶lonymi warto¶ciami. *)
    let minl = map selector
                 (filter 
		    (fun x -> value x = Undefined) 
		    (map (follow dag) e))
    in
      if e = [] && v = Undefined then 
	(* Wierzcho³ek nierozwiniêty. *) 
	1
      else if minl <> [] then 
	(* Zliczenie ilo¶ci wiercho³ków wyznaczaj±cych warto¶æ. *)
	cumulator minl
      else 
	(* Wszyscy synowie s± okre¶leni. *)
	infty

  (* Przelicz warto¶æ wierzcho³ka. *)
  and recalculate_value = 
    if v = Undefined && e <> [] then
      let ss = map (fun x -> value (follow dag x)) e
      in
        if who = And then 
	  fold_left and3 True ss
	else 
	  fold_left or3 False ss
    else v
  in
    (* Liczba wierzcho³ków falsyfikuj±cych. *)
    let minf = minft min_f (if who = Or then sum else minl)

    (* Liczba wierzcho³ków potwierdzaj±cych. *)
    and mint = minft min_t (if who = And then sum else minl)
    in
    
      update dag (s, who) (recalculate_value, minf, mint, e);;


(* Dodanie nierozwiniêtego li¶cia. *)
let insert_leaf dag s who = 
  if is_present dag s who then 
    dag 
  else 
    insert dag s who Undefined [];;

(* Pocz±tkowa sytuacja do rozwa¿enia. *)
let initial s = insert_leaf empty_map s i;;


(* Rozwija w danym dag-u wêze³ v odpowiadaj±cy sytuacji s dla gracza who. *)
let rec expand_dag dag s who = 

  (* Rozwija li¶æ v odpowiadaj±cy sytuacji s. *)
  let expand_leaf v = 
    if finished s who then begin
      (* Rozgrywka skoñczona. *)
      insert dag s who (result s who) []
    end else 
      (* Rozwiñ mo¿liwe ruchy. *)
      let leafs = map (fun m -> (m, move s who m)) (moves s who)
      in
        (* Wstaw nowych synów do DAG-u. *)
        let dag1 = 
	  fold_left 
	    (fun d (m, s) -> insert_leaf d s (other who)) 
	    dag leafs
        (* Lista nastêpników. *)
	and e = map (fun (m, s) -> (m, s, other who)) leafs
	in
	  insert dag1 s who (value v) e

  (* Przelicza warto¶æ wêz³a na podstawie warto¶ci potomków. *)
  and recalculate_value d v = 
    let ss = map (fun x -> value (follow d x)) (succ v)
    in
      let vv = 
	if who = And then 
	  fold_left and3 True ss 
	else 
	  fold_left or3 False ss
      in 
        insert d s who vv (succ v)
  in

    (* Rozwija podgraf zaczepiony w we¼le v wybieraj±c minimaln± *)
    (* warto¶æ podanego aspektu.                                 *)
    let expand_compound v aspect = 
      (* Predykat sprawdzaj±cy, czy syn wymaga rozwiniêcia. *)
      let test s = 
	let vv = follow dag s
	in aspect vv = aspect v && value vv = Undefined
      in
        (* Rozwijamy syna realizuj±cego min-f/min-t *)
        let ssl = filter test (succ v)
	in
	  if ssl = [] then 
	    (* Potomek jest ju¿ rozwiniêty. *)
	    recalculate_value dag v
	  else
	    let dag1 = expand_dag dag (succ_board (hd ssl)) (other who)
	    in
	      recalculate_value dag1 v
    in
      let v = vertex dag s who 
      in
        if not (is_expanded v) then 
	  (* Rozwiniêcie li¶cia. *)
	  expand_leaf v
	else if succ v = [] then 
	  (* Nie ma co rozwijaæ. *)
	  dag
	else if who = And then 
	  (* Rozwiniêcie wêz³a And. *)
	  expand_compound v min_f
	else
	  (* Rozwiniêcie wêz³a Or. *)
	  expand_compound v min_t;;


(* Rozwijaj, a¿ zostanie okre¶lona waga korzenia. *)
let rec keep_expanding dag s who = 
  if value (vertex dag s who) = Undefined then
    keep_expanding (expand_dag dag s who) s who
  else
    dag;;

(* Najlepszy ruch. *)
let game s = 
  let dag = initial s
  in
    let dag1 = keep_expanding dag s i
    in
      let v = vertex dag1 s i 
      in
        (succ_move (hd (filter 
			  (fun s -> value v = value (follow dag1 s))
			  (succ v))));;


