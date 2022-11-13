#load "dynlink.cma";;
#load "camlp4o.cma";;

open Lazy;;
open Camlp4.PreCast;;
open Syntax

(* Typ strumieniowy. *)
type 'a stream = Nil | Cons of 'a * 'a stream Lazy.t;;


(* Konstruktory i selektory strumieniowe *)
let empty s = s = Nil;;

let head s = 
  match s with
    Nil -> failwith "Empty" |
    Cons (h, _) -> h;;

let tail s = 
  match s with
    Nil -> failwith "Empty" |
    Cons (_, t) -> force t;; 

(* Makro Cons. *)
#load "pa_strum.cmo";;

#use "streams-impl.ml";;


