
let rotacja arr k = 
  let n = Array.length arr in
  let (+%) a b = (a + b) mod n in 
  let rem = arr.(0) in
  for i = 0 to n - 1 do begin
    arr.(i) <- arr.((n - k) +% i)
  end done; arr.(n - k) <- rem;;

rotacja [|'a', 'l', 'a', 'm', 'a', 'k', 'o', 't', 'a'|] 4;;