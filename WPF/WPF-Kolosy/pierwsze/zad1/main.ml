
(* KOLOKWIUM PIERWSZE - ZADANIE 1 *)
(* Aleksander TUDRUJ *)
(* 12.11.2020 *)


(********)
(* OPIS *)
(********)

(* funkcja podz li zwraca parę dwóch list. Pierwsza to lista elementów występujących wielokrotnie, druga to lista elementów występujących jednokrotnie *)
(* całe rozwiązanie polega na sprawdzaniu czy dla danego elementu listy, nazwijmy go x_i, x_i jest różne niż jego sąsiedzi (t.j. x_(i-1) oraz x_(x-2)) *)
(* rozpatruję wszystkie przypadki odpowiedinimi wzorcami *)
(* pierwszy element listy sprawdzany jest tylko z drugim, zaś ostatni tylko z przedostatnim *)
(* posortowanie elementów gwarantuje, że jezeli dla danego x_i ewentualne powtórzenia wartości będą na indeksach i-1 lub i+1 *)
(* spodziewana zlożoność czasowa to O(n) dla n=liczba_elemntów wejściowej listy *)
(* przechodzę się po tablicy dokładnie jeden raz, odwracanie wykonuję dwa razy, przy czym wykonam w sumie O(n) operacji *)

(******************)
(* OPIS ZMIENNYCH *)
(******************)

(* funkcja wysun -> lis -> lista podana na wejściu *)
(* funckja podz -> lista -> lista podana na wejściu *)
(* funkcja aux -> *)
(* li -> pozostała lista do rozpatrzenia *)
(* poprz -> poprzedni element *)
(* wielo_acc -> lista element występujących wielokrotnie (już znalezionych i stwierdzonych) *)
(* jedno_acc -> lista element występujących jednokrotnie (już znalezionych i stwierdzonych) *)

(*******)
(* KOD *)
(*******)

(* słowo: unikalne -> występujące dokładnie jeden raz w wejściowej liście *)

let wysun lis =
  let podz lista =
    let rec aux li poprz wielo_acc jedno_acc =
      match li with
      | [] -> (wielo_acc, jedno_acc) (* wszystko zostało już policzone *)
      | [x] -> if poprz = x then  (* osobno rozpatruję ostatni element z poprzednim *)
          (x :: wielo_acc, jedno_acc) (* x jest równe innemu elementowi *)
        else
          (wielo_acc, x :: jedno_acc) (* x jest unikalne w liście *)
      | h1 :: h2 :: t ->
        if h1 = poprz || h1 = h2 then (* czy x_i = x_(i-1) lub x_i = x_(i+1) *)
          aux (h2 :: t) h1 (h1 :: wielo_acc) jedno_acc
        else
          aux (h2 :: t) h1 wielo_acc (h1 :: jedno_acc)
    in match lista with
    | [] -> ([], [])
    | [x] -> ([], [x])
    | (h1 :: h2 :: t) -> if h1 = h2 then (* osobno rozpatruję głowę z następnym elementem *)
        aux (h2 :: t) h1 [h1] [] (* h1 nie jest unikalne *)
      else
        aux (h2 :: t) h1 [] [h1] (* h1 jest unikalne *)
  in let (wielo, jedno) = podz lis in
  (List.rev wielo) @ (List.rev jedno);; (* wrzucanie na akumulator odwróciło listy, więc odwracam je przed scaleniem *)

(*********)
(* TESTY *)
(*********)

wysun [1; 2; 2; 3; 4; 5; 5; 6; 6; 6; 7];;
wysun [1];;
wysun [];;
wysun [1; 2];;
wysun [1; 1];;
wysun [1; 1; 2];;
wysun [1; 2; 3; 4; 4; 4; 4];;
wysun [1; 1; 1; 2; 3; 4];;
wysun [1; 2; 2; 2; 4];;

wysun ["kot"; "pies"; "pies"; "zaba"; "zaba"; "zebra"];;
