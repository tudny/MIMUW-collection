
(* Graph with size and matrix of distances. *)
type graph = { size : int; distances : int array array }


let inf = 1000000009


let makeGraph n = 
  let dist = Array.make_matrix n n inf
  in 
  for i = 0 to n - 1 do
    Array.set (Array.get dist i) i 0
  done;
  { size = n;  distances = dist }


let getEdge { size = size; distances = dist } beg en = 
  Array.get (Array.get dist beg) en


let addEdge ({ size = size; distances = dist } as g) beg en len =
  let actualLen = getEdge g beg en in
  if actualLen > len then
    Array.set (Array.get dist beg) en len


let getN { size = size } = size


let readGraph () = 
  let n = ref 0 in
  let m = ref 0 in
  Scanf.scanf " %d %d" (fun x y -> n := x; m := y);
  let n = !n in
  let m = !m in
  let g = makeGraph n in

  for i = 0 to m - 1 do
    let a = ref 0 in
    let b = ref 0 in
    let l = ref 0 in
    Scanf.scanf " %d %d %d" (fun x y z -> a := x; b := y; l := z);
    let a = !a in
    let b = !b in
    let l = !l in

    assert (0 <= a && a < n);
    assert (0 <= b && b < n);
    assert (0 <= l && l < inf);

    addEdge g a b l
  done; g


let algorithmFloydWarshall g = 
  let n = getN g in
  let shortest = Array.make_matrix n n 0 in
  
  for i = 0 to n - 1 do
    for j = 0 to n - 1 do
      Array.set (Array.get shortest i) j (getEdge g i j)
    done
  done;

  let get arr x y = Array.get (Array.get arr x) y in

  for k = 0 to n - 1 do
    for i = 0 to n - 1 do
      for j = 0 to n - 1 do
        let newPath = get shortest i k + get shortest k j in
        if get shortest i j > newPath then
          Array.set (Array.get shortest i) j newPath
      done
    done
  done;

  shortest


let answerQuery shortest = 
  let q = ref 0 in
  Scanf.scanf " %d" (fun x -> q := x);
  let q = !q in 
  for i = 0 to q - 1 do
    let a = ref 0 in
    let b = ref 0 in
    Scanf.scanf " %d %d" (fun x y -> a := x; b := y);
    let a = !a in
    let b = !b in
    let anwser = Array.get (Array.get shortest a) b in
    if anwser == inf then
      Printf.printf "INF\n"
    else
      Printf.printf "%d\n" anwser
  done


let main () = 
  let g = readGraph () in
  let shortest = algorithmFloydWarshall g in
  answerQuery shortest
;;

main ()
;;