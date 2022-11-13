open List;;

(* Algorytm naiwny wyszukiwania wzorca. *)

let rec prefix w t = 
  match w with 
    [] -> true |
    hw::tw -> 
      match t with 
	[] -> false |
	ht::tt -> (hw = ht) && prefix tw tt;;
   
let pattern_match w t = 
  let rec iter t n acc = 
    let a = if prefix w t then n::acc else acc 
    in
      match t with 
	[] -> rev a |
	_::tail -> iter tail (n+1) a
  in 
    iter t 1 [];;

(* Algorytm KMP, wersja tablicowa. *)
open Array;;

(* Przygotowanie tablicy KMP *)
let pref t = 
  let
    p = make (length(t) + 1) 0 and
    pj = ref 0 
  in begin
    for i = 2 to length t do 
      while (!pj > 0) && (t.(!pj) <> t.(i - 1)) do
        pj := p.(!pj)
      done;
      if t.(!pj) = t.(i - 1) then pj := !pj + 1; 
      p.(i) <- !pj 
    done;
    p
  end;;

(* Zastosowanie KMP do wyszukiwania wzorca. *)
let find x y = 
  let 
    i = ref 0 and
    j = ref 0 and
    w = ref [] and
    p = pref x
  in 
    while !i <= length y - length x do
      j := p.(!j);
      while (!j < length x) && (x.(!j) = y.(!i + !j)) do 
        j := !j + 1
      done;
      if !j = length x  then w := !i::!w;
      i := !i + if !j > 0 then !j - p.(!j) else 1
    done;
    rev !w;;
    
let tekst = [| 'a'; 'l'; 'a'; 'l'; 'a'; 'b'; 'a'; 'm'; 'a' |];;
let w1 = [| 'a'; 'l'; 'a' |];;
let w2 = [| 'a' |];;
let w3 = [| 'x'; 'y'; 'z' |];;


(* Algorytm KMP, wersja kolejkowa. *)

#use "fifo-lifo.ml";;

let pref t = 
  let rec iter ti p j qj qi acc = 
    match ti with 
      [] -> rev acc |
      ht:: tt -> 
	if j > p then 
	  iter ti p (j - 1) (remove_last qj) (put_first qi (last qj)) acc
	else if (j > 0) && not (ht = snd (first qi)) then 
	  iter ti (fst (last qj)) j qj qi acc
	else if ht = snd (first qi) then 
	  iter tt (j+1) (j+1)  
	    (put_last qj (first qi)) 
	    (put_last (remove_first qi) (j+1, ht))
	    ((j+1)::acc)
	else 
	  iter tt 0 0 qj 
	    (put_last qi (0, ht)) 
	    (0::acc)
  in
    match t with 
      [] -> [] |
      ht::tt -> 
	iter tt 0 0 empty_queue (put_first empty_queue (0,ht)) [0];;

