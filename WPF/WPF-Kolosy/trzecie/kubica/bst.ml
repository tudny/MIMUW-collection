(* Zbiory skoñczone. *)
module type FINSET = 
  sig
    (* Zbiory skoñczone jako drzewa BST *)
    type 'a finset

    (* Wyj±tek podnoszony, gdy badamy warto¶æ spoza dziedziny. *)
    exception Undefined

    (* Pusty zbiór to puste drzewo *)
    val empty_finset : 'a finset

    (* Predykat charakterystyczny zbioru. *)
    val member : 'a finset -> 'a -> bool

    (* Wstawienie elementu do zbioru *)
    val insert : 'a finset -> 'a -> 'a finset
	
    (* Najmniejszy element *)
    val min : 'a finset -> 'a

    (* Najwiêkszy  element *)
    val max : 'a finset -> 'a

    (* Usuniêcie elementu ze zbiotu *)
    val remove : 'a finset -> 'a -> 'a finset

    (* Lista elementów zbioru. *)
    val finset2list : 'a finset -> 'a list
  end;;


module BST : FINSET = 
  struct
    (* Zbiory skoñczone jako drzewa BST *)
    type 'a finset = Empty | Node of ('a  * 'a finset * 'a finset)

    (* Wyj±tek podnoszony, gdy badamy warto¶æ spoza dziedziny. *)
    exception Undefined

    (* Pusty zbiór to puste drzewo *)
    let empty_finset = Empty

    (* Znajd¼ poddrzewo o zadanym kluczu              *)
    (* (je¶li nie ma, to wynikiem jest puste drzewo). *)
    let rec find s e =
      match s with
	Empty -> Empty |
	Node (k, l, r) ->
	  if e = k then s
	  else if e < k then find l e
	  else find r e

    (* Predykat charakterystyczny zbioru. *)
    let member s e = not (find s e = Empty)
    
    (* Wstawienie elementu do zbioru *)
    let rec insert s e =
      match s with
	Empty -> Node (e, Empty, Empty) |
	Node (k, l, r) ->
	  if e = k then
            s
	  else if e < k then
            Node (k, insert l e, r)
	  else
            Node (k, l, insert r e)

    (* Najmniejszy element *)
    let rec min s =
      match s with
	Node (e, Empty, _) -> e |
	Node (_, l, _)     -> min l |
	Empty              -> raise Undefined
    
    (* Najwiêkszy  element *)
    let rec max s =
      match s with
	Node (e, _, Empty) -> e |
	Node (_, _, r)     -> max r |
	Empty              -> raise Undefined
    
    (* Usuniêcie elementu z listy *)
    let rec remove s e = 
      match s with
	Empty -> s |
	Node (k, l, r) ->
	  if e = k then
            match (l, r) with
              (Empty, Empty) -> Empty |
              (Empty, _) ->
		let m = min r
		in Node (m, l, remove r m) |
		_  ->
		  let m = max l
		  in Node (m, remove l m, r)
	  else if e < k then
            Node (k, remove l e, r)
	  else
            Node (k, l, remove r e)

    (* Lista elementów zbioru. *)
    let finset2list s = 
      let res = ref [] 
      in
        let rec walk s = 
	  match s with 
	    Empty -> () |
	    Node (k, l, r) -> 
	      begin
		walk r;
		res := k :: !res; 
		walk l 
	      end
	in
          walk s;
          !res
  end;;
