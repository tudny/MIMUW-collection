
(* KOLOKWIUM PIERWSZE - ZADANIE 2 *)
(* Aleksander TUDRUJ *)
(* 12.11.2020 *)


(********)
(* OPIS *)
(********)

(* rozwiązanie polega na pamiętaniu trzech poprzednich elementów i na ich podstawie obliczanie nowego oraz wywołanie się rekurencyjne z nowymi elementami *)
(* pierwsze elementy (n \in {0, 1, 2}) są stałe *)
(* kolejne obliczane są ze wzoru *)
(* aby nie przechodzić za każdym razem ciągu b to pamiętam tylko aktulanie potrzebne elementy b (t.j. od jakiegoś indeksu do n-1) *)
(* reszta b wrzucana jest na akumulator i gdy dojdę do końca ciągu b odwracam akumulator i wstawiam znowu jako b (czyszcząc przy tym akumulator) *)
(* operacja odwracania ma złożoność czasową O(k) i wykonuję ją n/k razy (co k elementów muszę odnowić ciąg b) *)
(* zatem złożoność czasowa wynosi O(n/k * k) = O(n) (złożoność się amurtyzuje, jeżeli mogę tak powiedzieć) *)
(* na początku nie mogę mieć poczatku ciagu b, więc pierwsze 3 elementy od razy wrzucam na akumulator *)
(* pierwszy rozpatrywany element to taki o n=3, więc potrzebuję b_{3 mod k} (k>=3) co daje k_0 lub k_3 w zależności od k *)
(* tak naprawdę przechodzę się po ciagu b i robię to w kółko (dzięki akumulatorowi) *)
(* program zakończy się, pownieważ doiteruję się do szukanego elementu *)

(******************)
(* OPIS ZMIENNYCH *)
(******************)

(* bonifacy -> n, ciag_b z opisu zadania *)
(* aux *)
(* m3 -> jeżeli jesteśmy w elemencie i to m3 = a_{i-3} *)
(* m2 -> jeżeli jesteśmy w elemencie i to m2 = a_{i-2} *)
(* m1 -> jeżeli jesteśmy w elemencie i to m1 = a_{i-1} *)
(* szukany_n -> aktulanie poszukiwany element ciągu a (jego numer) *)
(* id -> aktualnie wyliczany element ciągu *)
(* ciag_b -> odpowiedni 'moment' przejścia po ciagu_b z treści zadania *)
(* acc_b -> już rozpatrzone elementy ciągu b, które później zostaną odwrócone i rozpatrywane ponowanie jako ciąg_b *)

(*******)
(* KOD *)
(*******)

let bonifacy n ciag_b =

  let rec aux m3 m2 m1 szukany_n id ciag_b acc_b =
    match ciag_b with
    | [] -> aux m3 m2 m1 szukany_n id (List.rev acc_b) [] (* jeżeli ciag_b jest pusty, to znaczy, że rozpatrzyliśmy już cały i możemy odwrócić akumulator *)
    | (b_imodk :: b_rest) ->
      if id = szukany_n + 1 then (* numer aktualnie rozpatrywanego elementu jest o jeden większy od tego, który mamy znaleźć, zatem m1 jest odpowiedzią *)
        m1
      else
        let akt_liczba = m1 + (if b_imodk = 0 then m2 else m3) in (* obliczam ze wzoru element a_{id} *)
        aux m2 m1 akt_liczba szukany_n (id + 1) b_rest (b_imodk :: acc_b) in (* dla następnego wywołania m2 stanie się m3, m1->m2, a_{id}(akt_liczba)->m1 *)

  let (h1 :: h2 :: h3 :: t) = ciag_b in (* ciąg_b z treści jest conajmniej trzyelementowy zatem mogę wyciągnąć z niego (bez błędu) trzy pierwsze elementy *)
  match n with
  | 0 -> 0
  | 1 -> 1
  | 2 -> 1
  | n -> aux 0 1 1 n 3 t [h3; h2; h1];; (* h1 to b_0, h2 to b_1, h3 to b_2, a dla id=3 potrzebne b_? to b_0 lub b_3 (w zależności od k) więc będzie ono albo pierwszym elementem t (pozostałości liczby b bez pierwszych trzech elementów) lub będzie to m1 po odwróceniu listy (co gwarantuje funkcja aux) *)
                                        (*aux (*a_0*)0 (*a_1*)1 (*a_2*)1 (*szukany_element*)n (*element_do_rozpatrzenia*)3 (*resztka ciągu_b*)t (*rozpatrzone elementy ciagu_b*)[h3; h2; h1];; *)

(*********)
(* TESTY *)
(*********)

bonifacy 0 [1; 1; 0; 1];;
bonifacy 1 [1; 1; 0; 1];;
bonifacy 2 [1; 1; 0; 1];;
bonifacy 3 [1; 1; 0; 1];;
bonifacy 4 [1; 1; 0; 1];;
bonifacy 5 [1; 1; 0; 1];;
bonifacy 6 [1; 1; 0; 1];;
bonifacy 7 [1; 1; 0; 1];;
bonifacy 8 [1; 1; 0; 1];;
bonifacy 9 [1; 1; 0; 1];;
bonifacy 10 [1; 1; 0; 1];;
bonifacy 11 [1; 1; 0; 1];;
bonifacy 12 [1; 1; 0; 1];;

bonifacy 8 [0; 1; 0; 1];;
bonifacy 9 [0; 1; 0; 1];;
