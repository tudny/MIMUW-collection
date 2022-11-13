let przeciecie (x1,x2) (y1,y2) =
  abs (min x2 y2 - max x1 y1);;

let powierzchnia ((x1,y1),(x2,y2)) ((x1p,y1p),(x2p,y2p)) = 
  przeciecie (x1,x2) (x1p,x2p) *
  przeciecie (y1,y2) (y1p,y2p);;
