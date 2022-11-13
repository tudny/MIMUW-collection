class punkt =
   object 
     val mutable x = 0
     val mutable y = 0
     method polozenie = (x,y)
     method przesun (xv, yv) = 
       begin
	 x <- x + xv;
	 y <- y + yv
       end
   end;;

let p = new punkt;;
p#przesun (1,0);;
p#przesun (3,2);;
p#polozenie;;

class punkt (x0, y0) =
   object 
     val mutable x = x0
     val mutable y = y0
     method polozenie = (x,y)
     method przesun (xv, yv) = 
       begin
	 x <- x + xv;
	 y <- y + yv
       end
     method na_poczatek = 
       begin
	 x <- x0;
	 y <- y0
       end
   end;;

let p = new punkt (8,5);;
p#przesun (-4, -3);;
p#polozenie;;
p#na_poczatek;;
p#polozenie;;

class kolorowy_punkt p c = 
  object
    val color : int = c
    inherit punkt p
    method kolor = color 
  end;;

let p = new kolorowy_punkt (3,1) 42;;
p#przesun (1,1);;
p#polozenie;;
p#kolor;;


(p : kolorowy_punkt :> punkt);;
(p :> punkt);;

class virtual punkt =
   object 
     val mutable x = 0
     val mutable y = 0
     method polozenie = (x,y)
     method virtual przesun : int * int -> unit
   end;;

class p = 
  object
    inherit punkt
    method przesun (xm, ym) = 
      begin
	x <- x + xm;
	y <- y + ym
      end
  end;;
