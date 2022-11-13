#load "camlp4o.cma";;
#load "pa_extend.cmo";;
#load "q_MLast.cmo";;

open Pcaml;;

(* Prosty uleniwiacz i poganiacz. *)
type 'a delayed = unit -> 'a;;


EXTEND            
  expr: LEVEL "simple"
  [ [ LIDENT "llazy"; param=expr -> 
    <:expr< let v () = $param$ in v >> ] ]
  ; 
END;;

let force x = x ();;

let a = llazy 4;;
force a;;

let duzo = llazy (print_string "6 x 9 = "; 42);;
force duzo;;


(* Efektywny uleniwiacz, poganiacz i opakowanie leniwej warto¶ci. *)
type 'a value = NoValue | Exception of exn | Value of 'a;;
let memoize f = 
  let 
    v = ref NoValue 
  in
    function () -> 
      match !v with 
        Value y -> y |
	Exception e -> raise e |
        NoValue -> 
          let 
            y = try f () 
	      with e -> begin
		v := Exception e;
		raise e
	      end
          in 
            v := Value y;
            y;;
           
let force x = x ();;


EXTEND            
  expr: LEVEL "simple"
  [ [ LIDENT "llazy"; param=expr -> 
    <:expr< let v () = $param$ in memoize v >> ] ]
  ; 
END;;

let a = llazy 4;;
force a;;

let b = llazy (print_int 42; 42);;
force b;;
force b;;



