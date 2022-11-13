(* GRAFY  - imperatywnie. *)


(* Modu� parametryzuj�cy ew. etykiety przypisywane kraw�dziom. *)
module type GRAPH_LABELS =
  sig
    type label
  end;;

(* Typ modu�u implementuj�cego grafy o wierzcho�kach numerowanych od 0 do n-1. *)
module type GRAPH =
  sig
    (* Etykiety kraw�dzi. *)
    type label

    (* Typ graf�w o wierzcho�kach typu 'a. *)
    type graph

    (* Pusty graf. *)
    val init : int -> graph

    (* Rozmiar grafu *)
    val size : graph -> int

    (* Dodanie kraw�dzi skierowanej ��cz�cej dwa wierzcho�ki. *)
    val insert_directed_edge : graph -> int -> int -> label -> unit

    (* Dodanie kraw�dzi nieskierowanej (dwukierunkowej) ��cz�cej dwa wierzcho�ki. *)
    val insert_edge : graph -> int -> int -> label -> unit

    (* Lista incydencji danego wierzcho�ka. *)
    val neighbours :  graph -> int -> (int*label) list
  end;;

(* Funktor tworz�cy grafy. *)
module Graph = functor (L : GRAPH_LABELS) ->
  (struct
    (* Etykiety kraw�dzi. *)
    type label = L.label

    (* Typ graf�w - tablica list s�siedztwa. *)
    type graph = {n : int; e : (int*label) list array}

    (* Pusty graf. *)
    let init s = {n = s; e = Array.make s []}

    (* Rozmiar grafu. *)
    let size g = g.n

    (* Dodanie kraw�dzi skierowanej ��cz�cej dwa wierzcho�ki. *)
    let insert_directed_edge g x y l =
      assert ((x<>y)  && (x >= 0) && (x < size g) && (y >= 0) && (y < size g) && (List.filter (fun (v,_) -> v=y) g.e.(x) = []));
      g.e.(x) <- (y,l)::g.e.(x)

    (* Dodanie kraw�dzi ��cz�cej dwa (istniej�ce) wierzcho�ki. *)
    let insert_edge g x y l =
      insert_directed_edge g x y l;
      insert_directed_edge g y x l

    (* Lista incydencji danego wierzcho�ka. *)
    let neighbours g x =
      g.e.(x)
  end : GRAPH with type label = L.label);;


#use "fifo.ml";;

(* Generyczny algorytm przeszukiwania. *)
module GraphSearch (Q : QUEUE) (G : GRAPH) =
  struct
    open G

    (* Obej�cie ca�ego grafu z wykonaniem procedury visit w wierzcho�kach, z przekazaniem jej innego koloru w ka�dej sk�adowej. *)
    let search visit g =
      let visited = Array.make (size g) false
      in
        (* Obej�cie sk�adowej grafu. *)
        let walk x color =
          let q = ref Q.empty
          in
            begin
              q := Q.insert !q x;
              visited.(x) <- true;
              while not (Q.is_empty !q) do
		let v = Q.front !q
		in begin
                  visit v color;
                  q := Q.remove !q;
                  List.iter
                    (fun (y, _) ->
                      if not (visited.(y)) then begin
			q := Q.insert !q y;
			visited.(y)<- true
                      end)
                    (neighbours g v);
		end
              done
            end
        in

          (* Kolor *)
          let k = ref 0
          in
            for i = 0 to size g - 1 do
              if not visited.(i) then begin
                (* Kolejna sk�adowa. *)
                walk i (!k);
                k := !k+1
              end;
            done;
            !k

  end;;


module DFS = GraphSearch (Lifo);;

module BFS = GraphSearch (Fifo);;

module FloatLabel =
  struct
    type label = float
  end;;

module type FLOAT_GRAPH = GRAPH with type label = float;;

module FloatGraph = Graph(FloatLabel);;

module Dijkstra (G : FLOAT_GRAPH) =
  struct
    open G

    module Order =
      struct

      end

    (* Obej�cie grafu z wierzcho�ka v, wynikiem jest tablica odleg�o�ci. *)
    let search g v =
      let visited = Array.make (size g) false
      and dist = Array.make (size g) infty
      and q = ref Q.empty
      in begin
        q := Q.insert !q x;
        visited.(x) <- true;
        while not (Q.is_empty !q) do
          visit (Q.front !q) color;
          List.iter
            (fun (y, _) ->
              if not (visited.(y)) then begin
                q := Q.insert !q y;
                visited.(y)<- true
              end)
            (neighbours g (Q.front !q));
          q := Q.remove !q
        done
      end
        in

             (* Kolor *)
          let k = ref 0
          in
            for i = 0 to size g do
              if not visited.(i) then begin
                (* Kolejna sk�adowa. *)
                walk i (!k);
                k := !k+1
              end;
            done;
            !k
  end;;












(**********************************************************************************)



(* GRAFY  - funkcyjnie. *)
#use "mapy.ml";;
#use "bst.ml";;


module Graph : GRAPH =
  struct
    (* Typ graf�w o wierzcho�kach typu 'a. *)
    type 'a graph = { v :  'a BST.finset ; e : ('a, 'a BST.finset) BstMap.map}

    (* Pusty graf. *)
    let empty = { v = BST.empty_finset; e = BstMap.empty}

    (* Dodanie wierzcho�ka (je�li istnieje, to b.z.). *)
    let insert_vertex g x =
      if BST.member g.v x then g else { v = BST.insert g.v x; e = BstMap.update g.e x BST.empty_finset }

    (* Dodanie kraw�dzi ��cz�cej dwa (istniej�ce) wierzcho�ki. *)
    let insert_edge g x y =
      assert (BST.member g.v x);
      assert (BST.member g.v y);
      assert (x <> y);
      let e1 = BstMap.update g.e x (BST.insert (BstMap.apply g.e x) y)
      in
        let e2 = BstMap.update e1 y (BST.insert (BstMap.apply e1 y) x)
        in { v = g.v; e = e2}


    (* Lista wszystkich wierzcho�k�w. *)
    let vertices g = BST.finset2list g.v

    (* Lista incydencji danego wierzcho�ka. *)
    let neighbours g x =
      BST.finset2list (BstMap.apply g.e x)

  end;;


(* Generyczny algorytm przeszukiwania. *)
module GraphSearch = functor (Q : QUEUE) ->
  struct
    open Graph

    let search visit g =
      let visited = ref BST.empty_finset
      in
        let walk x =
          let q = ref Q.empty
          in
            begin
              q := Q.insert !q x;
              visited := BST.insert !visited x;
              while not (Q.is_empty !q) do
                visit (Q.front !q);
                List.iter
                  (fun y ->
                    if not (BST.member !visited y) then begin
                      q := Q.insert !q y;
                      visited := BST.insert !visited y;
                    end)
                  (neighbours g (Q.front !q));
                q := Q.remove !q
              done
            end
        in
          List.iter (fun x -> if not (BST.member !visited x) then walk x) (vertices g)
  end;;

module DFS = GraphSearch (Lifo);;

module BFS = GraphSearch (Fifo);;

