(* Konto bankowe *)

(* Wersja bez enkapsulacji *)

let saldo = ref 0;;

let wplata kwota = 
  if kwota > 0 then 
    (saldo := !saldo + kwota;
    !saldo)
  else
    failwith "Z³a kwota";;
  
let wyplata kwota = 
  if kwota > 0 then 
    if kwota <= !saldo then 
      (saldo := !saldo - kwota;
      !saldo)  
    else 
      failwith "Brak ¶rodków na koncie"
  else 
    failwith "Z3a kwota";;
      

(* Wersja z enkapsulacj± *)

let zaloz poczatek = 
  let 
    saldo = ref poczatek
  in 
    let wplata kwota =
      if kwota > 0 then 
        (saldo := !saldo + kwota;
        !saldo) 
      else  
        failwith "Ujemna lub zerowa kwota"
    and wyplata kwota =
      if kwota > 0 then 
        if kwota <= !saldo then 
          (saldo := !saldo - kwota;
          !saldo) 
        else 
          failwith "Brak ¶rodksw na koncie"
      else 
        failwith "Ujemna lub zerowa kwota"
    in 
      (wplata, wyplata);;
      
