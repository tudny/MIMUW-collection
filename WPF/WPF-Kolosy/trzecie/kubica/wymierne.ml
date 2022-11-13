(* Liczby wymierne *)


(* Prosta implementacja struktury danych *)

let ulamek l m = (l, m);;

let licznik (l, m) = l;;

let mianownik (l, m) = m;;

(* Konstrukcja ulamka ze skracaniem, ale bez ujednolicenia znaków *)

let ulamek l m =
  let r = nwd l m
  in ((l / r), (m / r));;
 

(*========== B A R I E R A   A B S T R A K C J I ============*)

(* Modyfikatory i równo¶æ, implementacja niezale¿na od struktury danych *)

let add_wym x y = 
  ulamek 
    (licznik x * mianownik y + licznik y * mianownik x)
    (mianownik x * mianownik y);;

let sub_wym x y = 
  ulamek 
    (licznik x * mianownik y - licznik y * mianownik x)
    (mianownik x * mianownik y);;

let mul_wym x y = 
  ulamek 
    (licznik x * licznik y) 
    (mianownik x * mianownik y)

let div_wym x y = 
  ulamek 
    (licznik x * mianownik y) 
    (mianownik x * licznik y)

let equal_wym x y = 
  licznik x * mianownik y = licznik y * mianownik x;;

 