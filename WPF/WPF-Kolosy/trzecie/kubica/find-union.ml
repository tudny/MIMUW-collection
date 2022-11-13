(* Implementacja drzew find-union. *)
module type FIND_UNION = sig
  (* Typ klas elementów typu 'a. *)
  type 'a set

  (* Tworzy now± klasê z³o¿on± tylko z danego elementu. *)
  val make_set : 'a -> 'a set

  (* Znajduje reprezentanta danej klasy. *)
  val find : 'a set -> 'a

  (* Sprawdza, czy dwa elementy s± równowa¿ne. *)
  val equivalent : 'a set -> 'a set -> bool 

  (* Scala dwie dane (roz³±czne) klasy. *)
  val union : 'a set -> 'a set -> unit

  (* Lista elementów klasy. *)
  val elements : 'a set -> 'a list

  (* Liczba wszystkich klas. *)
  val n_of_sets : unit-> int
end;;


module Find_Union : FIND_UNION = struct
  (* Typ klas elementów typu 'a. *)
  type 'a set = { 
    elem         : 'a;           (* Element klasy. *)
    up           : 'a set ref;   (* Przodek w drzewie find-union. *)
    mutable rank : int;          (* Ranga w drzewie find-union. *)
    mutable next : 'a set list   (* Lista potomków w pewnym drzewie rozpinaj±cym klasê. *)
  }

  (* Licznik klas. *)
  let sets_counter = ref 0

  (* Liczba wszystkich klas. *)
  let n_of_sets () = !sets_counter

  (* Tworzy now± klasê z³o¿on± tylko z danego elementu. *)
  let make_set x = 
    let rec v = { elem = x; up = ref v; rank = 0; next = [] }
    in begin
      sets_counter := !sets_counter + 1;
      v
    end

  (* Znajduje korzeñ drzewa, kompresuj±c ¶cie¿kê. *)
  let rec go_up s = 
    if s == !(s.up) then s 
      else begin
	s.up := go_up !(s.up);
	!(s.up)
      end

  (* Znajduje reprezentanta danej klasy. *)
  let find s = 
    (go_up s).elem
	
  (* Sprawdza, czy dwa elementy s± równowa¿ne. *)
  let equivalent s1 s2 =
    go_up s1 == go_up s2

  (* Scala dwie dane (roz³±czne) klasy. *)
  let union x y = 
    let fx = go_up x 
    and fy = go_up y
    in
      if not (fx == fy) then begin
	if fy.rank > fx.rank then begin
	  fx.up := fy;
          fy.next <- fx :: fy.next
        end else begin
	  fy.up := fx;
	  fx.next <- fy :: fx.next;
	  if fx.rank = fy.rank then fx.rank <- fy.rank + 1
	end;
	sets_counter := !sets_counter - 1
      end
  
  (* Lista elementów klasy. *)
  let elements s = 
    let acc = ref []
    in
      let rec traverse s1 = 
	begin
	  acc := s1.elem :: !acc;
	  List.iter traverse s1.next
	end
      in begin
	traverse (go_up s);
	!acc
      end


  (* Testy *)
  let _ = 
    let n1 = n_of_sets()
    and s1 = make_set 42
    and s2 = make_set 24
    and s3 = make_set 12
    and s4 = make_set 21
    in  begin
      assert (n1 = 0);
      assert (n_of_sets() = 4);
      assert (find s1 <> find s2);
      assert (elements s1 = [42]);
      assert (elements s2 = [24]);
      assert (elements s3 = [12]);
      assert (elements s4 = [21]);
      union s1 s2;
      union s3 s4;
      assert (n_of_sets() = 2);
      assert (List.sort compare (elements s1) = [24; 42]);
      assert (List.sort compare (elements s2) = [24; 42]);
      assert (List.sort compare (elements s3) = [12; 21]);
      assert (List.sort compare (elements s4) = [12; 21]);
      union s4 s3; 
      assert (n_of_sets() = 2);
      assert (find s1 = find s2);
      assert (find s1 <> find s3);
      assert (find s3 = find s4);
      union s2 s4;
      assert (n_of_sets() = 1);
      assert (find s1 = find s3);
      assert (List.sort compare (elements s1) = [12; 21; 24; 42]);
      assert (List.sort compare (elements s2) = [12; 21; 24; 42]);
      assert (List.sort compare (elements s3) = [12; 21; 24; 42]);
      assert (List.sort compare (elements s4) = [12; 21; 24; 42]);
    end

end;;
