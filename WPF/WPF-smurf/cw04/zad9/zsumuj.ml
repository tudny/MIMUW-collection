
open List;;

let get_indeks (_, a, _) = a;;

let rec popraw_sume oczekiwany_indeks obiekt_sumy =
    if oczekiwany_indeks >= get_indeks obiekt_sumy then obiekt_sumy
    else
        let (suma, ost_indeks, (h :: t)) = obiekt_sumy in
        popraw_sume oczekiwany_indeks (suma + h, ost_indeks - 1, t);;


let zsumuj lista =
    let rec aux li id poprz_suma =
        if id = 0 then []
        else
            let (suma, ost_indeks, lista_left) = popraw_sume (hd li) poprz_suma in
            suma :: (aux (tl li) (id - 1) (suma, ost_indeks, lista_left))
    in rev (aux (rev lista) (length lista) (0, 1 + length lista, rev lista));;


let lista_do_sprawdzenia = [2; 4; 4; 5; 6; 6; 6; 10];;

zsumuj lista_do_sprawdzenia;;
zsumuj [1; 2; 3; 4];;
zsumuj [1; 2; 2; 2];;
zsumuj [2; 2; 2; 3; 6; 10];;
zsumuj [1; 2; 3; 4; 5; 6];;

assert (( zsumuj [1; 2; 3; 4] ) = [10; 9; 7; 4]);;
