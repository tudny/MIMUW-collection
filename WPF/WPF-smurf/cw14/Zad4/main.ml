(*
type pora = Wiosna | Lato | Jesien | Zima;;


let pory li =
  let rec worm p k score dis cnt0 cnt1 cnt2 cnt3 =
    print_int cnt0; print_string ", ";
    print_int cnt1; print_string ", ";
    print_int cnt2; print_string ", ";
    print_int cnt3; print_string " | ";
    print_int dis; print_string ", p:";
    print_int (List.length p); print_string ", k:";
    print_int (List.length k); print_string "\n";
    match p, k with
    | [], _ -> score
    | _, [] -> score
    | (ph :: pt, kh :: kt) ->
      if cnt0 > 0 && cnt1 > 0 && cnt2 > 0 && cnt3 > 0 then
        let f = worm pt k (min score dis) (dis - 1)  in
        match ph with
        | Wiosna -> f (cnt0 - 1) cnt1 cnt2 cnt3
        | Lato   -> f cnt0 (cnt1 - 1) cnt2 cnt3
        | Jesien -> f cnt0 cnt1 (cnt2 - 1) cnt3
        | Zima   -> f cnt0 cnt1 cnt2 (cnt3 - 1)
      else
        let f = worm p kt score (dis + 1) in
        if kt = [] then score else
        let next = List.hd kt in
        match next with
        | Wiosna -> f (cnt0 + 1) cnt1 cnt2 cnt3
        | Lato   -> f cnt0 (cnt1 + 1) cnt2 cnt3
        | Jesien -> f cnt0 cnt1 (cnt2 + 1) cnt3
        | Zima   -> f cnt0 cnt1 cnt2 (cnt3 + 1)
  in
  let wynik =
    let f = worm li li max_int 1 in
      match List.hd li with
      | Wiosna -> f 1 0 0 0
      | Lato   -> f 0 1 0 0
      | Jesien -> f 0 0 1 0
      | Zima   -> f 0 0 0 1
  in if wynik = max_int then -1 else wynik
;;

pory [Lato; Wiosna; Zima; Jesien; Lato; Zima; Wiosna; Jesien; Lato];; *)