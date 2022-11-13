
let heads family =
  List.flatten (List.map (fun li -> try [List.hd li] with | _ -> []) family)

let lista = [[1; 2; 3]; [4; 5; 6]; [1; 2]; [4; 5]; [4]; [17; 20]; []];;

heads lista;;
