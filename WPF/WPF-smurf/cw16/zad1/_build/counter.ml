
let global_timer = ref 0

type counter = int ref * int ref

let make () = 
  (ref 0, ref !global_timer)

let inc (v, time) = 
  if !time < !global_timer then begin
    v := 1; time := !global_timer
  end else begin
    v := !v + 1
  end; !v

let reset () = 
  global_timer := !global_timer + 1

