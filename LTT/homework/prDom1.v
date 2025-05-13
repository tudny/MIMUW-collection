(* Aleksander Tudruj - at429630 *)

(* ZADANIA NA ĆWICZENIA LUB DOMOWE *)

Section Simple_types.

(* Można je udowodnić, zamienić Admitted. na Qed. i obejrzeć termy
dowodowe *)

  Lemma imp_dist : (P -> Q -> R) -> (P -> Q) -> P -> R.
  Proof.
   intros.
   apply H.
   - apply H1.
   - apply H0.
     apply H1.
  Qed.

Print imp_dist.
(*
imp_dist =
fun (H : P -> Q -> R) (H0 : P -> Q) (H1 : P) => H H1 (H0 H1)
     : (P -> Q -> R) -> (P -> Q) -> P -> R
*)

  Lemma imp_trans : (P -> Q) -> (Q -> R) -> P -> R.
  Proof.
   intros.
   apply H0.
   apply H.
   exact H1.
  Qed.

  Print imp_trans.
(*
imp_trans = fun (H : P -> Q) (H0 : Q -> R) (H1 : P) => H0 (H H1)
     : (P -> Q) -> (Q -> R) -> P -> R
*)

  Lemma imp_perm : (P -> Q -> R) -> Q -> P -> R.
  Proof.
   intros.
   apply H.
   - exact H1.
   - exact H0.
  Qed.

  Print imp_perm.
(*
imp_perm = fun (H : P -> Q -> R) (H0 : Q) (H1 : P) => H H1 H0
     : (P -> Q -> R) -> Q -> P -> R
*)

  Lemma ignore_Q : (P -> R) -> P -> Q -> R.
  Proof.
   intros.
   apply H.
   exact H0.
  Qed.

  Print ignore_Q.
(*
ignore_Q = fun (H : P -> R) (H0 : P) (_ : Q) => H H0
     : (P -> R) -> P -> Q -> R
*)

  Lemma delta_imp : (P -> P -> Q) -> P -> Q.
  Proof.
   intros.
   apply H.
   - exact H0.
   - exact H0.
  Qed.

  Print delta_imp.
(*
delta_imp = fun (H : P -> P -> Q) (H0 : P) => H H0 H0
     : (P -> P -> Q) -> P -> Q
*)

  Lemma delta_impR : (P -> Q) -> P -> P -> Q.
  Proof.
   intros.
   apply H.
   assumption.
  Qed.

  Print delta_impR.
(*
delta_impR = fun (H : P -> Q) (_ H0 : P) => H H0
     : (P -> Q) -> P -> P -> Q
*)

  Lemma diamond : (P -> Q) -> (P -> R) -> (Q -> R -> S) -> P -> S.
  Proof.
   intros.
   apply H1.
   - apply H.
     exact H2.
   - apply H0.
     exact H2.
   Qed.

  Print diamond.
(*
diamond =
fun (H : P -> Q) (H0 : P -> R) (H1 : Q -> R -> S) (H2 : P) => H1 (H H2) (H0 H2)
     : (P -> Q) -> (P -> R) -> (Q -> R -> S) -> P -> S
*)

(* Typ bez inhabitanta w lambda rachunku z typami prostymi - jak to
sprawdzić *)
  Lemma peirce : (((P -> Q) -> P) -> P).
  Proof.
  Admitted.
(*
Really tried, but failed.
*)

  Lemma weak_peirce : ((((P -> Q) -> P) -> P) -> Q) -> Q.
   intros.
   apply H.
   intros.
   apply H0.
   intros.
   apply H.
   intros.
   exact H1.
  Qed.

Print weak_peirce.
(*
weak_peirce =
fun H : (((P -> Q) -> P) -> P) -> Q =>
H (fun H0 : (P -> Q) -> P => H0 (fun H1 : P => H (fun _ : (P -> Q) -> P => H1)))
     : ((((P -> Q) -> P) -> P) -> Q) -> Q
*)

End Simple_types.

