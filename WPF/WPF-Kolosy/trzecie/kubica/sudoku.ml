(* Program rozwi�zuj�cy �amig��wki sudoku. *)
#use "backtracking.ml";;

module Sudoku_Framework = 
  struct
    open List

    (* Tablicowa reprezentacja �amig��wki.              *)
    (* Tablica n^2 x n^2 wype�niona liczbami 0..n^2.    *)
    (* Puste pola s� reprezentowane przez zera.         *)
    type sudoku = int array array

    (* Rozmiar n danego sudoku. *)
    let size (s:sudoku) = truncate(sqrt (float (Array.length s)))

    (* Lista kolejnych liczb ca�kowitych od-do. *)
    let rec ints a b = if a > b then [] else a :: ints (a+1) b

    (* Przekszta�ca par� list w list� par -- odpowiednik produktu. *)
    let pairs l1 l2 = 
      flatten (map (fun i-> map (fun j -> (i,j)) l1) l2)

    (* Lista kolejnych indeks�w sudoku rozmiaru n. *)
    let indeksy s = 
      let n = size s 
      in ints 0 (n * n - 1)

    (* Lista cyfr, kt�rymi wype�niamy sudoku. *)
    let cyfry s = 
      let n = size s 
      in ints 1 (n * n)

    (* Lista par indeks�w reprezentuj�cych jeden "kwadrat". *)
    let kwadrat s =
      let n = size s
      in
        let s = ints 0 (n-1)
	in pairs s s

    (* Tworzy list� indeks�w do pustych miejsc w sudoku. *)
    let brakujace (s : sudoku) = 
      filter 
	(fun (i,j) -> s.(i).(j) = 0) 
	(pairs (indeksy s) (indeksy s))
  
    (* Sprawdza, czy cyfr� k mo�na wstawi� w miejscu (i,j). *)
    let valid_value (s : sudoku) i j k = 
      let valid_column i =  
	for_all (fun x -> s.(i).(x) <> k || x = j) (indeksy s)
      and valid_row j = 
	for_all (fun x -> s.(x).(j) <> k || x = i) (indeksy s)
      and valid_square i j =
	let n = size s 
	in
	  for_all
            (fun (x, y) -> (x, y) = (i,j) || s.(x).(y) <> k)
            (map (fun (x,y) -> (n*(i/n) + x, n*(j/n) + y)) (kwadrat s))
      in
        (valid_column i) && 
        (valid_row j) &&
        (valid_square i j)

    (* Wyznacza list� cyfr, kt�re mo�na wstawi� w miejscu (i,j). *)
    let valid_values (s : sudoku) i j = 
      filter (valid_value s i j) (cyfry s)

    (* Sprawdza, czy sudoku jest poprawnie wype�nione. *)
    let validate_sudoku s =
      (* Sortowanie cyfr. *)
      let sort_int l = 
	sort (fun x y -> if x>y then 1 else if x<y then -1 else 0) l
      in
        (* Sprawdza, czy i-ty wiersz jest permutacj� cyfr. *)
        let valid_row i = 
	  sort_int (map (fun j -> s.(i).(j)) (indeksy s)) = cyfry s
	and valid_column j = 
	  sort_int (map (fun i -> s.(i).(j)) (indeksy s)) = cyfry s
	and valid_square (i, j) =
	  let n = size s
	  in sort_int (map (fun (x,y) -> s.(n*i+x).(n*j+y)) (kwadrat s)) = cyfry s
	in
	  brakujace s = [] &&
	  for_all valid_row (indeksy s) &&
	  for_all valid_column (indeksy s) &&
	  for_all valid_square (kwadrat s)

    (* Wypisuje sudoku. *)
    let print_sudoku (s : sudoku) = 
      let szlaczek n = 
	begin
	  print_string "+";
	  for i = 1 to n do 
	    for j = 1 to n do 
	      print_string "---"
	    done;
	    print_string "-+"
	  done;
	  print_newline ()
	end
      and n = size s 
      in
        (* G�rna ramka. *)
        szlaczek n;
        (* Kolejne rz�dy kwadrat�w (z ramkami z bok�w i z do�u). *)
        for i = 0 to n-1 do 
	  (* Kolejne wiersze w rzedzie kwadrat�w. *)
	  for x = 0 to n-1 do
	    (* Ramka z lewej strony wiersza. *)
	    print_string"| ";
	    (* Kolejne kwadraty, przez kt�re przechodzi wiersz. *)
	    for j = 0 to n-1 do
	      (* Kom�rki. *)
	      for y = 0 to n-1 do 
		if  s.(n*i+x).(n*j+y) < 10 then print_string " ";
		print_int s.(n*i+x).(n*j+y);
		print_string " "
	      done;
	      (* Oddzielenie s�siaduj�cych kwadrat�w. *)
	      print_string "| "
	    done;
	    print_newline()
	  done;
	  (* Ramka oddzielaj�ca rz�dy kwadrat�w. *)
	  szlaczek n
	done;
      (* Odst�p w pionie. *)
      print_newline();
  end;;

(* Przyk�adowe �amig��wki do cel�w testowych. *)
let s0 : Sudoku_Framework.sudoku = 
  [| 
     [| 5; 6; 3; 1; 4; 2; 7; 9; 8 |];
     [| 1; 4; 2; 7; 9; 8; 3; 5; 6 |];
     [| 7; 8; 9; 6; 3; 5; 1; 2; 4 |];
     [| 9; 7; 4; 2; 8; 3; 6; 1; 5 |];
     [| 6; 2; 5; 9; 1; 7; 4; 8; 3 |];
     [| 8; 3; 1; 5; 6; 4; 2; 7; 9 |];
     [| 4; 5; 7; 3; 2; 9; 8; 6; 1 |];
     [| 3; 9; 6; 8; 7; 1; 5; 4; 2 |];
     [| 2; 1; 8; 4; 5; 6; 9; 3; 7 |]
   |];;


let s1 : Sudoku_Framework.sudoku = 
  [| 
     [| 4; 0; 0; 0; 0; 0; 6; 0; 5 |];
     [| 0; 0; 7; 9; 6; 0; 0; 1; 0 |];
     [| 0; 6; 0; 4; 0; 0; 8; 0; 7 |];
     [| 0; 8; 3; 0; 0; 6; 0; 0; 0 |];
     [| 0; 1; 0; 0; 0; 0; 0; 5; 0 |];
     [| 0; 0; 0; 2; 0; 0; 4; 6; 0 |];
     [| 7; 0; 6; 0; 0; 9; 0; 2; 0 |];
     [| 0; 5; 0; 0; 2; 4; 7; 0; 0 |];
     [| 8; 0; 9; 0; 0; 0; 0; 0; 6 |]
   |];;

let w1 : Sudoku_Framework.sudoku = 
  [|
    [|4; 9; 8; 7; 1; 2; 6; 3; 5|]; 
    [|5; 3; 7; 9; 6; 8; 2; 1; 4|];
    [|1; 6; 2; 4; 5; 3; 8; 9; 7|]; 
    [|2; 8; 3; 5; 4; 6; 9; 7; 1|];
    [|6; 1; 4; 8; 9; 7; 3; 5; 2|]; 
    [|9; 7; 5; 2; 3; 1; 4; 6; 8|];
    [|7; 4; 6; 1; 8; 9; 5; 2; 3|]; 
    [|3; 5; 1; 6; 2; 4; 7; 8; 9|];
    [|8; 2; 9; 3; 7; 5; 1; 4; 6|]
  |];;

let s2 : Sudoku_Framework.sudoku = 
  [| 
     [| 0; 0; 0; 0; 0; 1; 0; 0; 8 |];
     [| 0; 0; 4; 7; 8; 0; 0; 0; 0 |];
     [| 0; 1; 0; 0; 6; 0; 2; 0; 0 |];
     [| 0; 4; 0; 3; 0; 0; 0; 0; 0 |];
     [| 0; 2; 1; 0; 0; 0; 6; 7; 0 |];
     [| 0; 0; 0; 0; 0; 5; 0; 8; 0 |];
     [| 0; 0; 7; 0; 4; 0; 0; 9; 0 |];
     [| 0; 0; 0; 0; 3; 2; 1; 0; 0 |];
     [| 5; 0; 0; 9; 0; 0; 0; 0; 0 |]
   |];;

let w2 : Sudoku_Framework.sudoku =
  [|
    [|6; 9; 3; 2; 5; 1; 7; 4; 8|]; 
    [|2; 5; 4; 7; 8; 9; 3; 1; 6|];
    [|7; 1; 8; 4; 6; 3; 2; 5; 9|]; 
    [|8; 4; 5; 3; 7; 6; 9; 2; 1|];
    [|3; 2; 1; 8; 9; 4; 6; 7; 5|]; 
    [|9; 7; 6; 1; 2; 5; 4; 8; 3|];
    [|1; 3; 7; 6; 4; 8; 5; 9; 2|]; 
    [|4; 8; 9; 5; 3; 2; 1; 6; 7|];
    [|5; 6; 2; 9; 1; 7; 8; 3; 4|]
  |];;

let s3 : Sudoku_Framework.sudoku = 
  [| 
     [| 0; 0; 4; 5; 0; 0; 2; 0; 1 |];
     [| 0; 0; 8; 3; 0; 0; 0; 7; 0 |];
     [| 1; 0; 0; 0; 0; 0; 0; 0; 6 |];
     [| 0; 2; 0; 0; 0; 6; 0; 0; 0 |];
     [| 0; 0; 3; 0; 5; 0; 7; 0; 0 |];
     [| 0; 0; 0; 4; 0; 0; 0; 8; 0 |];
     [| 8; 0; 0; 0; 0; 0; 0; 0; 9 |];
     [| 0; 7; 0; 0; 0; 1; 6; 0; 0 |];
     [| 2; 0; 1; 0; 0; 4; 3; 0; 0 |]
   |];;

let w3 : Sudoku_Framework.sudoku = 
  [|
    [|7; 3; 4; 5; 6; 8; 2; 9; 1|]; 
    [|9; 6; 8; 3; 1; 2; 5; 7; 4|];
    [|1; 5; 2; 9; 4; 7; 8; 3; 6|]; 
    [|5; 2; 9; 8; 7; 6; 4; 1; 3|];
    [|4; 8; 3; 1; 5; 9; 7; 6; 2|]; 
    [|6; 1; 7; 4; 2; 3; 9; 8; 5|];
    [|8; 4; 6; 7; 3; 5; 1; 2; 9|]; 
    [|3; 7; 5; 2; 9; 1; 6; 4; 8|];
    [|2; 9; 1; 6; 8; 4; 3; 5; 7|]
  |];;

let smaxi : Sudoku_Framework.sudoku = 
  [| 
     [| 0; 0; 5; 0;11; 0; 0; 2; 4; 0; 0;14; 0;10; 0; 0 |];
     [| 0; 9; 0; 8; 0; 0; 6; 0; 0; 1; 0; 0; 2; 0; 7; 0 |];
     [|14; 0;11; 0; 0; 3; 0; 0; 0; 0; 6; 0; 0; 8; 0;15 |];
     [| 0; 1; 0; 0;13; 0; 0; 0;10; 0; 0;16; 0; 0;12; 0 |];
     [|11; 0; 0; 4; 0; 0; 5; 0; 0; 3; 0; 0;13; 0; 0; 6 |];
     [| 0; 0;12; 0; 0;15; 0; 0; 0; 0;13; 0; 0;16; 0; 0 |];
     [| 0;15; 0; 0; 7; 0; 0;11; 9; 0; 0;12; 0; 0; 1; 0 |];
     [| 6; 0; 0;13; 0; 0;12; 0; 0;11; 0; 0; 3; 0; 0;14 |];
     [|10; 0; 0; 3; 0; 0;14; 0; 0;16; 0; 0; 5; 0; 0; 8 |];
     [| 0; 2; 0; 0; 5; 0; 0; 8;13; 0; 0; 6; 0; 0; 9; 0 |];
     [| 0; 0;15; 0; 0;16; 0; 0; 0; 0;12; 0; 0; 6; 0; 0 |];
     [|13; 0; 0;11; 0; 0; 1; 0; 0;14; 0; 0; 7; 0; 0; 2 |];
     [| 0; 7; 0; 0;12; 0; 0; 4; 8; 0; 0; 9; 0; 0;15; 0 |];
     [| 8; 0;16; 0; 0;10; 0; 0; 0; 0; 3; 0; 0;11; 0; 5 |];
     [| 0;14; 0; 1; 0; 0; 9; 0; 0;15; 0; 0; 8; 0; 2; 0 |];
     [| 0; 0; 6; 0; 2; 0; 0;16;12; 0; 0; 1; 0; 3; 0; 0 |]
   |];;



(* Typ modu�u implementuj�cego wnioskowanie n.t. sudoku. *)
module type SUDOKU_REASONING = 
  sig
    (* Konfiguracja, to plansza i lista pustych p�l. *)
    type config =  Sudoku_Framework.sudoku * (int * int) list

    (* Wynikiem reason jest lista p�l wywnioskowanych i pozosta�ych pustych. *)
    val reason : config -> (int * int) list * (int * int) list
  end;;

(* Funktor, kt�ry wplata procedur� wnioskowania o sudoku *)
(* w backtracking znajduj�cy rozwi�zania.                *)
module Sudoku_BT (R : SUDOKU_REASONING) = 
  struct
    open Sudoku_Framework

    (* Instancja problemu dla back-trackingu. *)
    module Problem = 
      struct 
	(* Wynikiem jest plansza sudoku. *)
	type result = sudoku

        (* Konfiguracja, to plansza i lista pustych p�l. *)
	type config = R.config

        (* Sudoku jest sko�czone, gdy nie ma wolnych p�l. *)
	let final (s, l) = l = []

	(* Kopia planszy. *)
	let copy_sudoku s = 
	  Array.init (Array.length s) (fun i-> Array.copy s.(i))
	
	(* Extract tworzy kopi� planszy. *)
	let extract (s, l) = copy_sudoku s


    
	(* Procedura przegl�daj�ca konfiguracje osi�galne w jednym kroku. *)
	(* Wype�niane jest pierwsze pole z listy pustych p�l.             *)
	let iter p (s, l) = 
	  if not (final (s, l)) then 
	    let (wst, l1) = R.reason (s, l)
	    in begin
	      if not (final (s, l1)) then begin
		(* Wybierz z listy pozycji do wype�nienia t� o najmniejszej *)
		(* liczbie mo�liwych miejsc do wstawienia.                  *)
		let sorted_moves = 
		  let cmp (x, _, _) (y, _, _) = 
		    if x < y then -1 else if x > y then 1 else 0
		  in
		    List.sort cmp 
	              (List.map 
			 (fun (i, j) -> 
			   let v = valid_values s i j
			   in (List.length v, (i, j), v))
		         l1)
		in
	          (* Pole z najkrotsza list� mo�liwo�ci. *)
	          let (_, (i, j), v) = List.hd sorted_moves
	          (* Pozosta�e pola. *)
		  and t = List.map (fun (_, p, _) -> p) (List.tl sorted_moves)
		  in 
		    List.iter 
	              (fun k -> 
			  s.(i).(j) <- k; 
			  p (s, t); 
			  s.(i).(j) <- 0
			)
                      v
	      end else p (s, l1);
	      List.iter (fun (i,j) -> s.(i).(j) <- 0) wst
	    end
      end

    (* Instancja rozwi�zania za pomoc� back-trackingu. *)
    module Solution = Bt_Solver (Problem)

    (* Procedura znajduj�ca rozwi�zania dla danej planszy sudoku. *)
    let sudoku s = 
      Solution.solutions (s, brakujace s)

  end;;

(* Brak wnioskowania o sudoku. *)
module NoReasoning = 
  struct
    (* Konfiguracja, to plansza i lista pustych p�l. *)
    type config =  Sudoku_Framework.sudoku * (int * int) list

    (* Nic nie wnioskujemy. *)
    let reason (s,l)  = ([],l)
  end;;

module Sudoku1 = Sudoku_BT (NoReasoning);;

assert (Sudoku1.sudoku s0 = [s0]);;
assert (Sudoku1.sudoku s1 = [w1]);;
assert (Sudoku1.sudoku s2 = [w2]);;
assert (Sudoku1.sudoku s3 = [w3]);;
(* Dla smaxi dzia�a ju� zbyt d�ugo .... *)

module Reasoning = 
  struct
    open Sudoku_Framework

    (* Konfiguracja, to plansza i lista pustych p�l. *)
    type config =  sudoku * (int * int) list

    (* Sprawdza czy dana warto�� wyst�puje na li�cie. *)
    let present x l = List.filter (fun y -> x = y) l <> []

    (* Je�eli w wierszu jest tylko jedno miejsce, w kt�rym mo�e  *)
    (* wyst�powa� jaka� cyfra, to jest tam ona wstawiana.        *)
    let reason_row (s, l) = 
      let changed = ref true
      and puste = ref l
      and wstawione = ref []
      and n = size s
      in 
        while !changed do
	  changed := false;
	  (* Kolejne wiersze. *)
	  for i = 0 to n * n - 1 do 
	    (* Warto�ci, kt�re mog� pojawi� si� na poszczeg�lnych miejscach. *)
	    let v = 
	      Array.init (n*n) 
		(fun j -> if s.(i).(j) = 0 then valid_values s i j else [s.(i).(j)])
	    in 
              (* Przejrzyj kolejne cyfry. *)
              for k = 1 to n * n do 
		(* Sprawd�, na ilu pozycjach mo�e wyst�powa� cyfra. *)
		let jj = List.filter (fun j -> present k v.(j)) (indeksy s)
		in
		  (* Je�li pozycja cyfry jest wyznaczona *)
		  (* jednoznacznie, wype�nij j�.         *)
		  match jj with 
		    [j] -> if s.(i).(j) = 0 then begin
		      changed := true;
		      s.(i).(j) <- k;
		      wstawione := (i,j)::!wstawione
                    end |
		    _ -> ()
	      done
	  done;
	  (* Korekta listy pustych p�l. *)
	  puste := List.filter (fun (i,j) -> s.(i).(j)=0) (!puste);
	done;
        (!wstawione, !puste)


    (* Je�eli w kolumnie jest tylko jedno miejsce, w kt�rym mo�e  *)
    (* wyst�powa� jaka� cyfra, to jest tam ona wstawiana.         *)
    let reason_col (s, l) = 
      let changed = ref true
      and puste = ref l
      and wstawione = ref []
      and n = size s
      in 
        while !changed do
	  changed := false;
	  (* Prostsze wnioskowanie. *)
	  let (w, p) = reason_row (s, !puste)
	  in begin
	    wstawione := w @ !wstawione; 
	    puste := p
	  end;
	  (* Kolejne kolumny. *)
	  for i = 0 to n * n - 1 do 
	    (* Warto�ci, kt�re mog� pojawi� si� na poszczeg�lnych miejscach. *)
	    let v = 
	      Array.init (n*n) 
		(fun j -> if s.(j).(i) = 0 then valid_values s j i else [s.(j).(i)])
	    in 
              (* Przejrzyj kolejne cyfry. *)
              for k = 1 to n * n do 
		(* Sprawd�, na ilu pozycjach mo�e wyst�powa� cyfra. *)
		let jj = List.filter (fun j -> present k v.(j)) (indeksy s)
		in
		  (* Je�li pozycja cyfry jest wyznaczona *)
		  (* jednoznacznie, wype�nij j�.         *)
		  match jj with 
		    [j] -> if s.(j).(i) = 0 then begin
		      changed := true;
		      s.(j).(i) <- k;
		      wstawione := (j,i)::!wstawione
                    end |
		    _ -> ()
	      done
	  done;
	  (* Korekta listy pustych p�l. *)
	  puste := List.filter (fun (i,j) -> s.(i).(j)=0) (!puste);
	done;
        (!wstawione, !puste)


    (* Je�eli w malym kwadracie jest tylko jedno miejsce, w kt�rym         *)
    (* mo�e wyst�powa� jaka� cyfra, to jest tam ona wstawiana.             *)
    let reason (s,l) = 
      let changed = ref true
      and puste = ref l
      and wstawione = ref []
      and n = size s
      in 
        while !changed do
	  changed := false;
	  (* Prostsze wnioskowanie. *)
	  let (w, p) = reason_col (s, !puste)
	  in begin
	    wstawione := w @ !wstawione; 
	    puste := p
	  end;
	  (* Kolejne kwadraty. *)
	  for i = 0 to n - 1 do 
	    for j = 0 to n - 1 do 
            (* Warto�ci, kt�re mog� pojawi� si� na poszczeg�lnych miejscach. *)
            let v = 
	      Array.init n (fun x -> 
		Array.init n (fun y -> 
		  if s.(i*n+x).(j*n+y) = 0 then 
                    valid_values s (i*n+x) (j*n+y)
		  else 
                    [s.(i*n+x).(j*n+y)]))
            in 
              (* Przejrzyj kolejne cyfry. *)
              for k = 1 to n * n do 
		(* Sprawd�, na ilu pozycjach mo�e wyst�powa� cyfra. *)
		let xy = List.filter (fun (x, y) -> present k v.(x).(y)) (kwadrat s)
		in
		  (* Je�li pozycja cyfry jest wyznaczona *)
		  (* jednoznacznie, wype�nij j�.         *)
  	  	  match xy with 
                    [(x, y)] -> if s.(i*n+x).(j*n+y) = 0 then begin
                      changed := true;
                      s.(i*n+x).(j*n+y) <- k;
		      wstawione := (i*n+x,j*n+y)::!wstawione
                    end |
                _ -> ()
	      done
	    done   
	  done;
	  (* Korekta listy pustych p�l. *)
	  puste := List.filter (fun (i,j) -> s.(i).(j)=0) (!puste)
	done;
        (!wstawione, !puste)


  end;;

module Sudoku2 = Sudoku_BT (Reasoning);;

assert (Sudoku2.sudoku s0 = [s0]);;
assert (Sudoku2.sudoku s1 = [w1]);;
assert (Sudoku2.sudoku s2 = [w2]);;
assert (Sudoku2.sudoku s3 = [w3]);;
assert (List.length (Sudoku2.sudoku smaxi) = 17);;



    


