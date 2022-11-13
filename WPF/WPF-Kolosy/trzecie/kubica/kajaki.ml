open List;;
#use "fifo-lifo.ml";;
     
(* P R O B L E M    K A J A K O W Y  

   W kajaku mo¿e byæ jedna lub dwie osoby  
   Kajak = lista wag siedz±cych w nim osób

   Rozwi±zanie I
   Dla ka¿dego elementu maksymalnego dobieramy najwiêkszy, 
   który siê z nim mie¶ci.  *)
   
let kajaki l wyp = 
  (* Najgrubszego kajakarza nazwiemy grubasem. 
     Tych spo¶ród pozosta³ych kajakarzy, którzy siê mieszcz± 
     z nim w kajaku nazwiemy chudzielcami. 
     Pozosta³ych równie¿ nazwiemy grubasami.  *)

  (* Skoryguj podzia³ na chudych i grubych. *)
  let rec dobierz g ch = 
    if (is_empty_queue g) && (is_empty_queue ch) then (g, ch) else
    if is_empty_queue g then
      (make_queue [last ch], remove_last ch) else
    if queue_size g = 1 then (g, ch) else 
    if first g + last g <= wyp then 
      dobierz (remove_first g) (put_last ch (first g)) 
    else (g, ch)
  in
    (* Obsad¼ jeden kajak. *)
    let rec sadzaj gp chp acc =
      let (g, ch) = dobierz gp chp
      in
        if is_empty_queue g then acc else
	if is_empty_queue ch then 
	  sadzaj (remove_last g) ch ([last g]::acc)
	else
	  sadzaj 
	    (remove_last g) 
	    (remove_last ch) 
	    ([last g; last ch]::acc)
    in
      sadzaj (make_queue l) empty_queue [];;

kajaki [1; 2; 2; 3; 3; 3; 4; 6; 6; 8; 8; 8; 9; 9; 10] 10;;

    
(* Rozwi±zanie II
   Dla ka¿dego elementu minimalnego dobieramy najwiêkszy, 
   który siê z nim mie¶ci. *)
   
let kajaki l wyp = 
  (* Najchudszy kajakarz jest chudzielcem. 
     Grubasy, to ci, którzy nie mieszcz± siê z nim w kajaku. 
     Pozostali to chudzielce. 
     Je¶li jest tylko jeden chudy, to przyjmujemy, ¿e jest on gruby. *)

  (* Skoryguj podzia³ na chudych i grubych. *)
  let rec dobierz ch g = 
    if is_empty_queue ch then (ch, g) else 
    if queue_size ch = 1 then 
      (empty_queue, put_first g (last ch)) else
    if first ch + last ch > wyp then 
      dobierz (remove_last ch) (put_first g (last ch))
    else 
      (ch, g)
in
  (* Obsad¼ jeden kajak. *)
  let rec sadzaj chp gp acc = 
    let (ch, g) = dobierz chp gp 
    in
      if (is_empty_queue ch) && (is_empty_queue g) then acc else
      if is_empty_queue ch then 
        sadzaj ch (remove_first g) ([first g]::acc)
      else 
        sadzaj (remove_first (remove_last ch)) g
               ([first ch; last ch]::acc)
  in
    sadzaj (make_queue l) empty_queue [];;

kajaki [1; 2; 2; 3; 3; 3; 4; 6; 6; 8; 8; 8; 9; 9; 10] 10;;


(* Rozwi±zanie III
   Je¶li siê da, to ³±czymy najgrubszego z najchudszym. *)
   
let kajaki l wyp = 
  let rec sadzaj q acc = 
    if is_empty_queue q then acc else 
    if queue_size q = 1 then [first q]::acc else 
    if first q + last q <= wyp then 
      sadzaj (remove_first (remove_last q)) ([first q; last q]::acc) 
    else 
      sadzaj (remove_last q) ([last q]::acc)
  in
    sadzaj (make_queue l) [];;

kajaki [1; 2; 2; 3; 3; 3; 4; 6; 6; 8; 8; 8; 9; 9; 10] 10;;


    
(* Rozwi±zanie IV
   Maksymalne dwa kolejne *)
   
let kajaki l wyp = 
  let rec dobierz ch g = 
    match (ch, g) with
      (_, []) -> (ch, []) |
      ([], h::t) -> dobierz [h] t |
      (chh::cht, gh::gt) -> 
        if chh + gh <= wyp then 
          dobierz (gh::ch) gt
        else 
          (ch, g)
  in 
    let rec sadzaj chp gp acc = 
      let (ch, g) = dobierz chp gp
      in 
        match ch with 
          [] -> acc |
          [h] -> sadzaj [] g ([h]::acc) |
          h1::h2::t -> sadzaj t g ([h2; h1]::acc)
    in 
      sadzaj [] l [];;

kajaki [1; 2; 2; 3; 3; 3; 4; 6; 6; 8; 8; 8; 9; 9; 10] 10;;



























