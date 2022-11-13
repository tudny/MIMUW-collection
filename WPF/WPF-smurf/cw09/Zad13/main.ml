
let prextrema lista =
  match lista with
  | [] -> []
  | (h :: t) ->
    let (_, _, ac) =
    List.fold_left
      (fun (mini, maks, acc) ele ->
         let _acc =
           if ele < mini || ele > maks then ele :: acc
           else acc
          in (min mini ele, max maks ele, _acc))
      (h, h, [h]) t
    in List.rev ac;;

prextrema [-2; 1; 0; 1; 3; 2; -1; 5; 4; -3; 2; 1; 7];;
