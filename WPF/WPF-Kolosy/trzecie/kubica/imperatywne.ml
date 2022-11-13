(* Programowanie imperatywne*)

(* Generator obiektów ze stanami lokalnymi *)
let generator s = 
  let r = ref s 
  in 
    function x -> begin
      r := !r + x;
      !r
    end;;
let o = generator 0;;
o 22;;
o 20;;

(* Konto bankowe *)
let konto pocz = 
  let saldo = ref pocz
  in
    let wplata kwota = 
      if kwota > 0 then begin 
	saldo := !saldo + kwota;
	!saldo
      end else failwith "Ujemna kwota"
    and wyplata kwota = 
      if kwota > 0 then 
	if kwota <= !saldo then begin 
	  saldo := !saldo - kwota;
	  !saldo
	end else failwith "Brak ¶rodków na koncie"
      else failwith "Ujemna kwota"
    in
      (wplata, wyplata);;

let (wplac, wyplac) = konto 0;;

wplac 50;;
wyplac 8;;

(* Zadanie zagadka: okre¶l co jest wynikiem tej procedury i udowodnij to. *)
let coto n = 
  let x = ref 0
  and y = ref 1
  in begin
    while !y <= n do 
      x := !x - !y;
      y := !y - !x
    done;
    if n mod 2 = 1 then !x else  (!x)
  end;;
    
