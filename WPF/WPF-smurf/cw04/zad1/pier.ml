
let pier n = 
    let rec aux id acc = 
        if id = 0 then acc
        else aux (id - 1) (id :: acc)
    in aux n [];;
    
pier 100;;
