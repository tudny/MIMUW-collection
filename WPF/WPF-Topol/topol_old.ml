
(* *********************************************** *)
(*                     WARNING                     *)
(*             THIS CODE DEOSN'T WORK              *)
(* *********************************************** *)


(* ******************************* *)
(* Topol (Sortowanie topologiczne) *)
(* Autor: Aleksander Tudruj        *)
(* Code review:                    *)
(* ******************************* *)

(* ***************** *)
(* Szkic rozwiązania *)
(* ***************** *)

(* Dla każadego wierzchołka zliczam ile krawędzi do niego wchodzi   *)
(* Tę liczbę nazwę stopniem wejścia                                 *)
(* Jeżeli dla danego wierzchołka stopień wejścia wynosi 0, to mogę  *)
(* dodać go do listy wynikowej oraz zmniejszyć stopień wejścia jego *)
(* dzieci, co rozumieć możemy jako usunięcie go z grafu             *)
(* Jeżeli stopień wejścia jakiegoś dziecka spadł do 0 to dodaję je  *)
(* do kolejki rozpatrywanych wierzchołków                           *)
(* Można zaoberwować, że jeżeli w grafie wystąpi cykl,              *)
(* to wierzchołki w tym cyklu nie zostaną nigdy dodane do kolejki   *)
(* skuktuje to niedodaniem ich do listy wynikowej                   *)
(* Zatem nie muszę robić testu na obecność cyklu w grafie;          *)
(* wystarczy sprawdzić czy lista wynikowa ma tyle wierzchołków      *)
(* ile jest ich w grafie                                            *)
(* Nie muszę również sprawdzać spójności grafu, gdyż nie wybieram   *)
(* Żadnego wierzchołka startowego, lecz konsekwentnie przechodzę    *)
(* się na początku po każdnym z nich                                *)

exception Cykliczne
;;

(* Funkcja przyjmuje graf w postaci danej na wejściu w zadaniu   *)
(* Zwraca mapę, której klucze i wartości odpowiadają numerowi    *)
(* wierzchołka oraz jego wartości podanej na wejściu             *)
(* Najprościej mówiąc mapa pozwala na ponumerowanie wierzchołków *)
let create_map graph = 
  let (nodes, neighbours) = List.split graph in
  let list_of_nodes = nodes @ List.flatten neighbours in
  List.fold_left
    (fun (_acc, cnt) node ->
      if PMap.mem node _acc then
        (_acc, cnt)
      else 
        (PMap.add node cnt _acc, cnt + 1)
    ) (PMap.empty, 0) list_of_nodes

(* Funkcja przyjmuje graf w postaci danej na wejściu w zadaniu    *)
(* oraz mapę wierzchołków (nazwaną wyżej numeracją)               *)
(* Zwraca graf w postaci list, lecz z przenazwanymi wierzchołkami *)
let convert graph mapped =
  List.map (fun (label, neighbours) ->
      (PMap.find label mapped,
       List.map (fun neighbour ->
          PMap.find neighbour mapped) neighbours)
    ) graph
;;

let prepare nr_of_nodes arr_graph = 
  let proper_graph = Array.init nr_of_nodes (fun i -> (i, [])) in
  Array.iter
    (fun (node, neighbours) -> 
      let (_, _actual) = proper_graph.(node) in
      proper_graph.(node) <- (node, neighbours @ _actual)
    ) arr_graph;
  proper_graph
;;

(* Funkcja przyjmuje graf w postaci danej na wejściu w zadaniu *)
(* oraz tablicę, która dla każdego wierzchołka będzie pamiętać *)
(* ile wierzchołków do niego wchodzi                           *)
let fill_in_degree graph in_degree =
  Array.iter (fun (label, neighbours) ->
      List.iter (fun neighbour ->
          in_degree.(neighbour) <- in_degree.(neighbour) + 1
        ) neighbours
    ) graph
;;

(* Funkcja przyjemuje kolejkę wierzchołków do przerobienia *)
(* oraz tablicę, która dla każdego wierzchołka pamięta ile *)
(* wierzchołków do niego wchodzi                           *)
let fill_queue queue in_degree =
  Array.iteri
    (fun id degree ->
      if degree = 0 then
        Queue.add id queue)
    in_degree
;;

(* Funkcja przyjmuje listę dzieci wierzchołka X i dla każdego dziecka *)
(* zmniejsza jego stopień wejścia (tablica [in_degree])               *)
let dec_degree neighbours in_degree =
  List.iter
    (fun neighbour -> in_degree.(neighbour) <- in_degree.(neighbour) - 1)
    neighbours
;;

(* Funkcja dodaje do kolejki rozpatrywanych wierzchołków wszystkie *)
(* wierzchołki, których stopień wejścia wynosi 0                   *)
let make_up_queue neighbours in_degree queue =
  List.iter
    (fun neighbour ->
       if in_degree.(neighbour) = 0 then Queue.add neighbour queue)
    neighbours
;;

(* Funkcja przyjmuje mapę i zwraca nową mapę powstałą przez zebranie *)
(* wszystkich par (key, value) i dodanie ich jako (value, key)       *)
let reverse_map map = 
  PMap.foldi
      (fun key value _acc -> PMap.add value key _acc)
      map PMap.empty
;;

(* Funkcja ta na podstawie wyniku, liczby wierzchołków oraz    *)
(* zmapowanego grafu tworzy wynik programu                     *)
(* Sprawdza czy wynik jest poprawny. Jeżeli wynik jest         *)
(* niepoprawny, to podnosi wyjątek [Cykliczne]                 *)
(* Jeżeli wynik jest poprawny przemapowuje numery wierzchołków *)
(* na ich oryginalnego wartości (etykietki)                    *)
let produce_output score map nr_of_nodes = 
  if List.length score < nr_of_nodes then
    raise Cykliczne
  else
    List.fold_left
      (fun _acc el -> PMap.find el map :: _acc) 
      [] score
;;

let get_neighbours arr_graph node = 
  arr_graph.(node) |> snd
;;

(* Funcja główna wykonująca listę kroków rozwiązania *)
let topol graph =
  (* Stworenie ponumerowania wierzchołków *)
  let (mapped, nr_of_nodes) = create_map graph in
  (* Przekonwertowanie wierzchołków z etykietek na numery *)
  let converted = convert graph mapped in
  (* Konwersja Listy wierzchołków na Tablicę *)
  (* Pozwala to na szybki dostęp do konkretnego wierzchołka *)
  let arr_graph = Array.of_list converted |> prepare nr_of_nodes in
  (* Tablica 'stopnia wejścia' dla każdego wierzchołka           *)
  (* 'Stopniem wejścia' danego wierzchołka nazwę liczbę krawędzi *)
  (* wchodzących do danego wierzchołka                           *)
  (* Z początku pusta; zaraz potem wypełniona przez przejście    *)
  (* się po wszystkich wierzchołkach i ich dzieciach             *)
  let in_degree = Array.make nr_of_nodes 0 in
  fill_in_degree arr_graph in_degree;
  (* Kolejka wierzchołków do rozpatrzenia                 *)
  (* Przez rozpatrzenie rozumiem usinięcie go z grafu     *)
  (* co skutkuje zmiejszeniem stopnia wejścia jego dzieci *)
  let queue = Queue.create () in
  fill_queue queue in_degree;
  (* Wynik zapisuję na liście (dodając do niej imperatywnie) *)
  let score = ref [] in
  let add_to_score x = score := x :: !score in
  (* Na liście znajdują się wierzchołki, których stopień wejścia wynosi 0 *)
  (* zatem mogę usunąć je z grafu, zmniejszyć stopień wejścia ich synów   *)
  (* oraz dodać na listę wynikową                                         *)
  while not (Queue.is_empty queue) do
    let element = Queue.take queue in
    let neighbours = get_neighbours arr_graph element in
    add_to_score element;
    dec_degree neighbours in_degree;
    make_up_queue neighbours in_degree queue
  done;
  (* Wynik przetwarzam oraz zwracam *)
  let reversed_map = reverse_map mapped in
  produce_output !score reversed_map nr_of_nodes
;;


