(* krok lista i max_ele duża_lista mała_lista*)

let rev_every lista =
  let rec rev_e li acc =
    match li with
    | [] 	 -> acc
    | h :: t -> rev_e t ((List.rev h) :: acc)
  in rev_e lista [];;


let podziel lista =
  let rec krok li id maks_ele akt_acc big_acc =
    match li with
    | [] 	 -> big_acc
    | h :: t ->
      let akt_maks = max h maks_ele in
      if id = akt_maks then krok t (id + 1) akt_maks [] ((h :: akt_acc) :: big_acc)
      else krok t (id + 1) akt_maks (h :: akt_acc) big_acc
  in rev_every (krok lista 1 0 [] []);;

podziel [1; 2; 3];;
podziel [1; 3; 2];;
podziel [2; 3; 1];;
podziel [1; 2; 3; 4; 5; 6; 7; 8];;
podziel [9; 8; 7; 6; 5; 4; 3; 2; 1];;
podziel [3; 2; 4; 6; 1; 5; 8; 7];;

podziel [1; 2; 3; 4; 5; 9; 8; 6; 7; 10; 11; 13; 12];;
