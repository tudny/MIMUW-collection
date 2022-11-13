(* Kodowanie pary int-ów na int-ie. *)
 let para a b =
   let rec pom a b c p =
     if (a = 0) && (b = 0) then c
     else pom (a/10) (b/10) (c + p *(a mod 10) + 10 * p *(b mod 10)) (100*p)
   in pom a b 0 1;;
 
 let rzut_x x =
   let rec pom x p a =
     if x = 0 then a
     else pom (x / 100) (10 * p) (a + (x mod 10) * p)
   in pom x 1 0;;
 
 let rzut_y x =
   let rec pom x p a =
     if x = 0 then a
     else pom (x / 100) (10 * p) (a + ((x/10) mod 10) * p)
   in pom x 1 0;;
 
