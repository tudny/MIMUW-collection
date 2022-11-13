
 type expr = 
  | And of expr * expr
  | Or of expr * expr
  | Not of expr
  | Value of bool;;


(* dp (true, false) *)
let rec wartosciowania e =
  let rec loop = function 
    | Value v -> (1, 1)
    | And (v, w) ->
      let (tl, fl) = loop v
      and (tr, fr) = loop w in
      (tl * tr, fl * tr + fr * tl + fl * fr)
    | Or (v, w) -> 
      let (tl, fl) = loop v
      and (tr, fr) = loop w in
      (tl * fr + tr * fl + tr * tl, fl * fr)
    | Not w ->
      let (t, f) = loop w in (f, t)
  in fst (loop e)
;;

wartosciowania (Or (Value true, Value false));;
wartosciowania (And((Or (Value true, Value false)), (Not(And(Or (Value false, Value false), And(Value false,Value false))))));;
wartosciowania (Or ((And((Or (Value true, Value false)), (Not(And(Or (Value false, Value false), And(Value false,Value false)))))), (And((Or (Value true, Value false)), (Not(And(Or (Value false, Value false), And(Value false,Value false))))))))

