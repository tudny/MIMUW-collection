(* Kolejki dwustronne ze sklejaniem i odwracaniem. *)
(* Implementacja imperatywna.                      *)

module type QUEUE = 
  sig
    (* Typ kolejki. *)
    type 'a queue

    (* Wyj±tek podnoszony, gdy kolejnka jest pusta, a nie powinna. *)
    exception EmptyQueue

    (* Konstruktor pustej kolejki. *)
    val init : unit -> 'a queue

    (* Predykat sprawdzaj±cy czy kolejka jest pusta. *)
    val is_empty : 'a queue -> bool

    (* Wstawienie na pocz±tek. *)
    val put_first : 'a queue -> 'a -> unit

    (* Wstawienie na koniec. *)
    val put_last : 'a queue -> 'a -> unit

    (* Pierwszy element. *)
    val first : 'a queue -> 'a

    (* Ostatni element. *)
    val last : 'a queue -> 'a

    (* Usuniêcie pierwszego el. *)
    val remove_first : 'a queue -> unit

    (* Usuniêcie ostatniego el. *)
    val remove_last : 'a queue -> unit

    (* Sklejenie dwóch kolejek, druga kolejka ulega wyczyszczeniu. *)
    val merge : 'a queue -> 'a queue -> unit

    (* Odwrócenie kolejki. *)
    val rev : 'a queue -> unit
  end;;

    
module Queue : QUEUE = 
  struct
    (* Typ kolejki. *)
    (* Kolejka sk³ada siê z rekordów u³o¿onych w listê dwukierunkow±,   *)
    (* ze stra¿nikami na koñcach; przy czym nie jest okre¶lone które    *)
    (* z dowi±zañ l1 i l2 prowadzi w przód, a które w ty³.              *)
    (* Wszystkie elementy poza stra¿nikami zawieraj± okre¶lon± warto¶æ. *)
    (* W przypadku stra¿ników dowi±zanie l1 wskazuje na stra¿nika.      *)
    (* Kolejka to para dowi±zañ do stra¿ników na koñcach.               *)
    type 'a opt = Null | Val of 'a
    type 'a elem = { mutable l1 : 'a elem; mutable l2 : 'a elem; v : 'a opt }
    type 'a queue = { mutable front : 'a elem; mutable back : 'a elem }

    (* Wyj±tek podnoszony, gdy kolejnka jest pusta, a nie powinna. *)
    exception EmptyQueue

    (* Warto¶æ elementu -- nie mo¿e byæ stra¿nik. *)
    let value e = 
      match e.v with 
	Val x -> x |
	Null  -> raise EmptyQueue

    (* Predykat sprawdzaj±cy czy kolejka jest pusta. *)
    let is_empty q = 
      (* Kolejka jest pusta gdy zawiera tylko stra¿ników. *)
      let f = q.front 
      and b = q.back 
      in
        assert (not (f == b));
        assert (f.l1 == f);
        assert (b.l1 == b);
        assert ((f.l2 == b) = (b.l2 == f));
        (f.l2 == b)

    (* Konstruktor pustej kolejki. *)
    let init () = 
      let rec g1 = { l1 = g1; l2 = g2; v = Null}
      and     g2 = { l1 = g2; l2 = g1; v = Null}
      in 
        { front = g1; back = g2 }

    (* Wstawia nowy element z warto¶ci± x pomiêdzy elementy p i q. *)
    let put_between p q x = 
      let r = { l1 = p; l2 = q; v = Val x }
      in begin
	assert (not (p == q) && 
		((p.l1 == q) || (p.l2 == q)) &&
	        ((q.l1 == p) || (q.l2 == p)));
	if p.l1 == q then p.l1 <- r else p.l2 <- r;
	if q.l1 == p then q.l1 <- r else q.l2 <- r
      end

    (* Wstawienie na pocz±tek. *)
    let put_first q x =
      let f = q.front 
      in
        put_between f f.l2 x

    (* Wstawienie na koniec. *)
    let put_last q x =
      let b = q.back
      in
        put_between b b.l2 x

    (* Pierwszy element. *)
    let first q = 
      value q.front.l2
	    
    (* Ostatni element. *)
    let last q = 
      value q.back.l2

    (* Usuwa element (nie bêd±cy stra¿nikiem). *)
    let remove e = 
      assert (not (e.l1 == e) && (e.v <> Null));
      let e1 = e.l1 
      and e2 = e.l2 
      in begin
	assert (not (e1 == e2) && not (e1 == e) && not (e2 == e));
        if e1.l1 == e then e1.l1 <- e2 else e1.l2 <- e2;
        if e2.l1 == e then e2.l1 <- e1 else e2.l2 <- e1
      end

    (* Usuniêcie pierwszego el. *)
    let remove_first q = 
      if is_empty q then 
	raise EmptyQueue
      else
	remove q.front.l2

    (* Usuniêcie ostatniego el. *)
    let remove_last q = 
      if is_empty q then 
	raise EmptyQueue
      else
	remove q.back.l2

    (* Sklejenie dwóch kolejek, druga kolejka ulega destrukcji. *)
    let merge q1 q2 = 
      assert (not (q1 == q2));
      let f1 = q1.front
      and b1 = q1.back
      and f2 = q2.front
      and b2 = q2.back
      in
        assert (not (f1==b1) && not (f2==b2));
        let e1 = b1.l2
	and e2 = f2.l2
	in begin
	  (* Z³±czenie kolejek. *)
	  if e1.l1 == b1 then e1.l1 <- e2 else e1.l2 <- e2;
	  if e2.l1 == f2 then e2.l1 <- e1 else e2.l2 <- e1;
	  (* Wynik w q1 *)
	  q1.back <- b2;
	  (* q2 puste *)
	  (* f2.l1 <- f2; *)
	  f2.l2 <- b1;
	  (* b1.l1 <- b1; *)
	  b1.l2 <- f2;
	  q2.back <- b1
	end

    (* Odwrócenie kolejki. *)
    let rev q = 
      let f = q.front
      and b = q.back 
      in begin
	q.front <- b;
	q.back <- f
      end

    (* Testy *)
    let _ = 
      begin
	assert (is_empty (init ()));
      end

  end;;

    
