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

(* Tr�jwarto�ciowa logika: prawda, fa�sz i warto�� nieokre�lona. *)
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

(* Gra w kamienie - cz�� zale�na od gry *)
#use "gra-w-kamienie.ml";;

    
(* Drzewo AND-OR jest reprezentowane za pomoc� mapy przyporz�dkowuj�cej *)
(* parom (sytuacja, gracz) n-tki postaci:                               *)
(* - warto��: True / False / Undefined                                  *)
(* - minimalna liczba el. falsyfikuj�cych,                              *)
(* - minimalna liczba el. potwierdzaj�cych,                             *)
(* - dla rozwini�tych: lista tr�jek (ruch, sytuacja, gracz).            *)


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

(* W�ze� jest rozwini�ty, je�eli ma syn�w lub ma warto��. *)
let is_expanded v = (value v <> Undefined) || (succ v <> []);;


(* Operacje na DAG-u *)

(* Czy ju� jest taki wierzcho�ek? *)
let is_present dag s who = dom dag (s, who);;

(* Dodanie do DAG-u w�z�a o parametrach:         *)
(* - s   - sytuacja,                             *)
(* - who - gracz,                                *)
(* - v   - wartosc,                              *)
(* - e   - lista tr�jek (ruch, sytuacja, gracz). *)
let insert dag s who v e =

  (* Rozmiar najmniejszego zbioru falsyfikuj�cego/potwierdzaj�cego *)
  (* (w zale�no�ci od parametr�w selector/cumulator) dla syn�w.    *)
  let minft selector cumulator =
    (*  Lista liczb do przeliczenia, pochodz�cych od  *)
    (* syn�w z nieokre�lonymi warto�ciami. *)
    let minl = map selector
                 (filter 
		    (fun x -> value x = Undefined) 
		    (map (follow dag) e))
    in
      if e = [] && v = Undefined then 
	(* Wierzcho�ek nierozwini�ty. *) 
	1
      else if minl <> [] then 
	(* Zliczenie ilo�ci wiercho�k�w wyznaczaj�cych warto��. *)
	cumulator minl
      else 
	(* Wszyscy synowie s� okre�leni. *)
	infty

  (* Przelicz warto�� wierzcho�ka. *)
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
    (* Liczba wierzcho�k�w falsyfikuj�cych. *)
    let minf = minft min_f (if who = Or then sum else minl)

    (* Liczba wierzcho�k�w potwierdzaj�cych. *)
    and mint = minft min_t (if who = And then sum else minl)
    in
    
      update dag (s, who) (recalculate_value, minf, mint, e);;


(* Dodanie nierozwini�tego li�cia. *)
let insert_leaf dag s who = 
  if is_present dag s who then 
    dag 
  else 
    insert dag s who Undefined [];;

(* Pocz�tkowa sytuacja do rozwa�enia. *)
let initial s = insert_leaf empty_map s i;;


(* Rozwija w danym dag-u w�ze� v odpowiadaj�cy sytuacji s dla gracza who. *)
let rec expand_dag dag s who = 

  (* Rozwija li�� v odpowiadaj�cy sytuacji s. *)
  let expand_leaf v = 
    if finished s who then begin
      (* Rozgrywka sko�czona. *)
      insert dag s who (result s who) []
    end else 
      (* Rozwi� mo�liwe ruchy. *)
      let leafs = map (fun m -> (m, move s who m)) (moves s who)
      in
        (* Wstaw nowych syn�w do DAG-u. *)
        let dag1 = 
	  fold_left 
	    (fun d (m, s) -> insert_leaf d s (other who)) 
	    dag leafs
        (* Lista nast�pnik�w. *)
	and e = map (fun (m, s) -> (m, s, other who)) leafs
	in
	  insert dag1 s who (value v) e

  (* Przelicza warto�� w�z�a na podstawie warto�ci potomk�w. *)
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

    (* Rozwija podgraf zaczepiony w we�le v wybieraj�c minimaln� *)
    (* warto�� podanego aspektu.                                 *)
    let expand_compound v aspect = 
      (* Predykat sprawdzaj�cy, czy syn wymaga rozwini�cia. *)
      let test s = 
	let vv = follow dag s
	in aspect vv = aspect v && value vv = Undefined
      in
        (* Rozwijamy syna realizuj�cego min-f/min-t *)
        let ssl = filter test (succ v)
	in
	  if ssl = [] then 
	    (* Potomek jest ju� rozwini�ty. *)
	    recalculate_value dag v
	  else
	    let dag1 = expand_dag dag (succ_board (hd ssl)) (other who)
	    in
	      recalculate_value dag1 v
    in
      let v = vertex dag s who 
      in
        if not (is_expanded v) then 
	  (* Rozwini�cie li�cia. *)
	  expand_leaf v
	else if succ v = [] then 
	  (* Nie ma co rozwija�. *)
	  dag
	else if who = And then 
	  (* Rozwini�cie w�z�a And. *)
	  expand_compound v min_f
	else
	  (* Rozwini�cie w�z�a Or. *)
	  expand_compound v min_t;;


(* Rozwijaj, a� zostanie okre�lona waga korzenia. *)
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


