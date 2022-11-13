
open Counter;;

let c1 = make ();;

let v1 = inc c1;;
Printf.printf "c1: %d\n" v1;;
let v2 = inc c1;;
Printf.printf "c1: %d\n" v2;;

let c2 = make ();;

let v3 = inc c2;;
Printf.printf "c2: %d\n" v3;;

reset ();;

Printf.printf "reset\n"

let v4 = inc c1;;
Printf.printf "c1: %d\n" v4;;

let v5 = inc c2;;
Printf.printf "c2: %d\n" v5;;