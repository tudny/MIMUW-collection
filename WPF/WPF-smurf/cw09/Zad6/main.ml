
let ostatki family =
  let last_singleton li = try [List.hd (List.rev li)] with | _ -> [] in
  List.flatten (List.map (last_singleton) family)

let lista = [[1; 2]; [1; 2; 3; 4]; [17; 20]; []; [40; 50]];;

ostatki lista;;
