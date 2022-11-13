
let tails l =
    let rec aux akt acc = 
        match akt with
        | [] -> [] :: acc
        | h :: t -> aux t (akt :: acc)
    in aux l [];;


tails [1; 2; 3];;

