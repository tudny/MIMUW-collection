
let cykl arr =
  let maks = ref 0 in
  let l = Array.length arr in
  let odw = Array.make l false in
  let rec dfs v acc =
    if odw.(v) = true then acc
    else begin
      odw.(v) <- true;
      dfs (arr.(v)) (acc + 1)
    end in
  for i = 0 to l - 1 do
    if odw.(i) = false then maks := max !maks (dfs arr.(i) 0)
  done; !maks;;


  cykl [|2; 1; 0; 5; 6; 4; 3; 8; 7|]



