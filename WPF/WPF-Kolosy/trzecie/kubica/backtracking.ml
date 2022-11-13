(********************************************************)
(* Ró¿ne generyczne sposoby rozwi±zania back-trackingu. *)
(********************************************************)

(* Generyczny opis problemu, który bêdziemy rozwi±zywaæ *)
(* za pomoc± back-trackingu.                            *)
module type BT_PROBLEM = 
  sig
    (* Typ rozpatrywanych konfiguracji. *)
    type config

    (* Typ szukanych wyników. *)
    type result
     
    (* Czy dana konfiguracja jest ju¿ kompletnym rozwi±zaniem. *)
    val final : config -> bool

    (* Wyci±ga z konfiguracji istotne rozwi±zanie.            *)
    (* W przypadu imperatywnych struktur danych tworzy kopiê. *)
    val extract : config -> result

    (* Procedura, która wywo³uje dan± procedurê, dla      *)
    (* dla ka¿dej konfiguracji osi±galnej w jednym kroku. *)
    val iter : (config -> unit) -> config -> unit
  end;;

(* Mechanizm szukaj±cy pierwszego lub wszystkich rozwi±zañ. *)
module type BT_SOLVER = functor (Problem : BT_PROBLEM) ->
  sig
    (* Wyj±tek podnoszony, gdy brak rozwi±zañ. *)
    exception NoSolution

    (* Procedura znajduj±ca jedno rozwi±zanie. *)
    val onesolution : Problem.config -> Problem.result

    (* Procedura znajduj±ca wszystkie rozwi±zania. *)
    val solutions : Problem.config -> Problem.result list
  end;;

(* Implementacja *)
module Bt_Solver : BT_SOLVER = functor (Problem : BT_PROBLEM) -> 
  struct
    (* Wyj±tek podnoszony, gdy brak rozwi±zañ. *)
    exception NoSolution

    (* Wyj±tek podnoszony, gdy znaleziono rozwi±zanie. *)
    exception Solution of Problem.result

    (* Procedura znajduj±ca rozwi±zanie. *)
    let onesolution s = 
      let rec backtrack s = 
	begin
	  if Problem.final s then 
	    raise (Solution (Problem.extract s));
	  Problem.iter backtrack s
	end
      in 
        try (backtrack s; raise NoSolution) 
	with Solution x -> x

    (* Procedura znajduj±ca wszystkie rozwi±zania. *)
    let solutions s = 
      let wynik = ref []
      in
        let rec backtrack s = 
	  begin
	    if Problem.final s then
	      wynik := (Problem.extract s)::!wynik;
	    Problem.iter backtrack s
	  end
	in begin      
	  backtrack s;
	 !wynik 
	end
  end;;

(********************************************************************)

(* Generyczny opis problemu, który bêdziemy rozwi±zywaæ *)
(* za pomoc± back-trackingu. Wersja z optymalizacj±.    *)
module type OBT_PROBLEM = 
  sig
    include BT_PROBLEM
     
    (* Dolne oszacowanie kosztu rozwi±zania,     *) 
    (* jakie mo¿na uzyskaæ z danej konfiguracji. *)
    val min_cost : config -> float

    (* Koszt kompletnego rozwi±zania. Musi byæ dodatni. *)
    val cost : config -> float
  end;;

(* Mechanizm szukaj±cy najtañszego rozwi±zania. *)
module type OBT_SOLVER = functor (Problem : OBT_PROBLEM) ->
  sig
    (* Wyj±tek podnoszony, gdy brak rozwi±zañ. *)
    exception NoSolution

    (* Procedura znajduj±ca rozwi±zanie. *)
    val solver : Problem.config -> Problem.result
  end;;

module OBt_Solver : OBT_SOLVER = functor (Problem : OBT_PROBLEM) -> 
  struct
    (* Wyj±tek podnoszony, gdy brak rozwi±zañ. *)
    exception NoSolution

    (* Opcjonalna warto¶æ. *)
    type 'a value = NoValue | Value of 'a

    (* Procedura znajduj±ca rozwi±zanie. *)
    let solver s = 
      let wynik = ref NoValue
      and koszt = ref infinity
      in
        let rec backtrack s = 
	  let c = Problem.cost s 
	  in begin
	    if (Problem.final s) && (c < !koszt) then begin
	      wynik := Value (Problem.extract s);
	      koszt := c
	    end;
	    if Problem.min_cost s < !koszt then 
	      Problem.iter backtrack s
	  end
	in begin      
	  backtrack s;
	  match !wynik with
	    Value x -> x |
	    NoValue -> raise NoSolution
	end
  end;;


