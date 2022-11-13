(* Testowanie strumieni. *)
let sample s = 
  let rec pom i s = 
    if i = 0 then []
    else head s :: pom (i-1) (tail s)
  in 
    pom 42 s;;


(* Dodanie sta³ej do strumienia *)
let rec add n s = 
  if empty s then Nil else Cons ((head s) + n, add n (tail s));;


(* Strumien kolejnych liczb ca³kowitych od n w górê. *)
let rec s n = 
  Cons (n, begin print_int n; print_string "\n"; s (n+1) end);;

let i = s 1;;
head (tail i);;
head (tail i);;


(* Strumieñ sta³y *)
let rec const_stream c = 
  Cons (c, const_stream c);; 

  
(* Strumieñ jedynek *)
let ones = const_stream 1;;


(* Filtrowanie strumienia *)
let rec filter_stream p s = 
  if empty s then 
    Nil
  else 
    let h = head s 
    and t = tail s
    in 
      if p h then 
        Cons (h, filter_stream p t)
      else
        filter_stream p t;;


(* Wyd³ubanie n-tego elementu strumienia *)
let rec stream_ref s n = 
  if n = 0 then 
    head s
  else
    stream_ref (tail s) (n - 1);;


(* Mapa dla strumieni *)    
let rec stream_map f s = 
  if empty s then 
    Nil
  else
    Cons (f (head s), stream_map f (tail s));;

(* Odpowiednik fold_left *)
let rec stream_fold f a s = 
  if empty s then 
    Cons (a, Nil)
  else
    Cons (a, stream_fold f (f a (head s)) (tail s));;

(* Odpowiednik fold-a dla strumieni -- rozwin±æ w przysz³o¶ci *)
let rec fold_stream f g a s = 
  if empty s then 
    a 
  else
    g a (head s) (lazy (fold_stream f g (f a (head s)) (tail s)));;

(* Filtrowanie strumieni z foldem *)
let filter_stream p s = 
  fold_stream 
    (fun a _ -> a)
    (fun _ h tl ->
      if p h then 
	Cons (h, (force tl))
      else
	force tl)
    Nil 
    s;;

(* Iterator *)    
let rec for_each p s = 
  if not (empty s) then begin
    p (head s);
    for_each p (tail s)
  end;;
  
let print_int_stream s = 
  let p x = 
    begin 
      print_int x;
      print_string "\n";
      flush stdout
    end
  in 
    for_each p s;;

let print_float_stream s = 
  let p x = 
    begin 
      print_float x;
      print_string "\n";
      flush stdout
    end
  in 
    for_each p s;;

    
(* Strumieñ ogonów. *)
let rec tails s = 
  Cons (s, tails (tail s));;    


(* Strumieñ liczb pierwszych. *)
let rec integers_from x = 
  Cons(x, integers_from (x+1));;

let nats = integers_from 0;;
  
let divisible x y = x mod y = 0;;

let sito p s = 
  filter_stream (function x -> not (divisible x p)) s;;
  
let rec sieve s  = 
  Cons (head s, sieve (sito (head s) (tail s)));;
  
let primes = 
  sieve (integers_from 2);;

print_int_stream primes;;


(* Uwik³ane definicje strumieni. *)
let rec ones =  
  Cons (1, ones);;

let rec add_int_streams s1 s2 = 
  if empty s1 or empty s2 then 
    Nil
  else
    Cons (head s1 + head s2, add_int_streams (tail s1) (tail s2));;

let rec integers =  
  Cons (1, add_int_streams ones integers);;

(* Strumieñ ogonów -- definicja uwik³ana. *)
let tails s = 
  let rec w = 
    Cons (s, stream_map tail w)
  in w;;

(* Liczby Fibonacciego - nie uwik³ane *)
let fibs = 
  let rec fibo a b = 
    Cons (a, fibo b (a+b))
  in fibo 0 1;;

(* Liczby Fibonacciego - uwik³ane *)
let rec fibs = 
  Cons (0, Cons (1, add_int_streams  fibs (tail fibs)));;
  
(* Uwik³ana definicja liczb pierwszych. *)
let square x = x * x;;

let rec primes = 
   Cons (2, filter_stream prime (integers_from 3))
and prime n = 
  let rec iter ps = 
    if square (head ps) > n then 
      true
    else if divisible n (head ps) then false
    else iter (tail ps)
  in 
    iter primes;;
    
    
(* Iteracje jako strumienie *)

let sqrt_improve x g = (g +. x /. g) /. 2.0;;

let sqrt_stream x = 
  let rec guesses = 
    Cons (1.0, stream_map (sqrt_improve x) guesses)
  in guesses;;
      
print_float_stream (sqrt_stream 2.0);;


(* Przybli¿enia pi *)

let pi_summands = 
  let 
    succ x = 
      if x > 0.0 then 
        -.x -. 2.0
      else
        -.x +. 2.0
  in
    let rec s = 
      Cons (1.0, stream_map succ s)
    in
      stream_map (fun x -> 1.0 /. x) s;;

let rec add_float_streams s1 s2 = 
  if empty s1 or empty s2 then 
    Nil
  else
    Cons (head s1 +. head s2, add_float_streams (tail s1) (tail s2));;
  
let partial_sums s = 
  let rec ps = 
      Cons (head s, add_float_streams (tail s) ps)
  in ps;;
   
let scale_stream c s = 
  stream_map (function x -> x *. c) s;;
  
let pi_stream = scale_stream 4.0 (partial_sums pi_summands);;
     
let square x = x *. x;;

let euler_transform st =
  let transform s = 
    let s0 = stream_ref s 0
    and s1 = stream_ref s 1
    and s2 = stream_ref s 2
    in (s2 -. square (s2 -. s1) /. (s0 -. 2.0 *. s1 +. s2))
  in
    stream_map transform (tails st);;

print_float_stream (euler_transform pi_stream);;

let rec pi_table = 
  Cons (pi_stream, stream_map euler_transform pi_table);;    

let pi_stream_acc = stream_map head pi_table;;

let pi = stream_ref pi_stream_acc 8;;



(* Scalenie dwóch monotonicznych strumieni liczb *)
let rec merge s1 s2 =
  if empty s1 then s2 
  else if empty s2 then s1 
  else 
    let h1 = head s1 
    and h2 = head s2
    in 
      if h1 < h2 then Cons (h1, merge (tail s1) s2)
      else if h1 > h2 then Cons (h2, merge s1 (tail s2))
      else Cons (h1, merge (tail s1) (tail s2));;


(* Scalenie dwóch strumieni uporz±dkowanych zgodnie z funkcj± wagi *)
let rec merge_weighted w s1 s2 = 
  if empty s1 then s2 else
  if empty s2 then s1 else 
    let h1 = head s1 
    and h2 = head s2
    in 
      if w h1 < w h2 then Cons (h1, merge_weighted w (tail s1) s2)
      else if w h1 > w h2 then Cons (h2, merge_weighted w s1 (tail s2))
      else Cons (h1, Cons (h2, merge_weighted w (tail s1) (tail s2)));;


(* Przekszta³ca strumieñ elementów uporz±dkowanych wg. wag w *)
(* strumieñ wag, które powtarzaj± siê                        *)
let rec non_uniq w s = 
  let rec skip we s =
    if empty s then Nil 
    else if we = w (head s) then skip we (tail s)
    else s
  in
    if empty s then Nil 
    else if empty (tail s) then Nil 
    else
      let h1 = head s
      and h2 = head (tail s)
      in 
        if w h1 = w h2 then 
	  Cons (w h1, non_uniq w (skip (w h1) s))
	else
	  non_uniq w (tail s);;

(* Strumieñ par (nieuporz±dkowanych) liczb naturalnych uporz±dkowany *)
(* zgodnie z sumami sze¶cianów wspó³rzêdnych                         *)
let ramanujan = 
  let cube x = x * x * x
  in
    let weight (x, y) = cube x + cube y
    in 
      let rec pairs s = 
	Cons ((head s, head s), 
              merge_weighted 
		weight
		(stream_map (function x -> (head s, x)) (tail s))
		(pairs (tail s)))
      in 
        non_uniq weight (pairs nats);;




(* Trójk±t Pascala *)
let rec pascal = 
  let next s = 
    let rec wyn = 
      Cons (1, add_int_streams (tail s) wyn)
    in 
      wyn
  in 
    Cons (ones, stream_map next pascal);;

List.map sample (sample pascal);;


(* Nieskoñczone s³owo Fibonacciego. *)

(* Rozwi±zanie oparte na obserwacji, ¿e a * Fib_1 * Fib_2 * .... * Fib_k = Fib_{k+2}. *)
let fibos = 
  let rec gen s = 
    if head s = 'a' then Cons ('a', Cons ('b', gen (tail s)))
    else Cons ('a', gen (tail s))
  in 
    let rec s = Cons ('b', gen s)
    in Cons ('a', s);;
  
(* Rozwi±zanie przez sklejanie kolejnych s³ów Fibonacciego. *)
let rec lapp s1 ls2 = if empty s1 then Lazy.force ls2 else Cons (head s1, lapp (tail s1) ls2);;

let fibos2 = 
  let rec fibs n =
    if n=1 then Cons ('b', Nil)
    else if n=2 then Cons ('a', Nil)
    else lapp (fibs (n-1)) (lazy (fibs (n-2)))
  and fibsy n = lapp (fibs n) (lazy (fibsy (n+1)))
  in lapp (fibs 2) (lazy (fibsy 1));;

assert (sample fibos = sample fibos2);;
