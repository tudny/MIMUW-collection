Section propositional_logic.

Variables P Q R S : Prop.

Lemma and_assoc : P /\ (Q /\ R) -> (P /\ Q) /\ R.
Proof.
  intro.
  destruct H.
  destruct H0.
  split.
  - split; assumption.
  - assumption.
Qed.

Lemma and_imp_dist : (P -> Q) /\ (R -> S) -> P /\ R -> Q /\ S.
Proof.
  intros.
  destruct H.
  destruct H0.
  split.
  - apply H.
    assumption.
  - apply H1.
    assumption.
Qed.


Lemma de_morgan_1 : ~(P \/ Q) -> ~P /\ ~Q.
Proof.
  unfold not.
  intros.
  split.
  - intro.
    apply H.
    left.
    assumption.
  - intro.
    apply H.
    right.
    assumption.
Qed.

Lemma de_morgan_2 : ~P /\ ~Q -> ~(P \/ Q).
Proof.
  unfold not.
  intros.
  destruct H.
  destruct H0.
  - apply H.
    assumption.
  - apply H1.
    assumption.
Qed.

Lemma de_morgan_3 : ~P \/ ~Q -> ~(P /\ Q).
Proof.
  unfold not.
  intros.
  destruct H0.
  destruct H.
  - apply H.
    assumption.
  - apply H.
    assumption.
Qed.

(* de_morgan_4 is not a tautology of IPL *)
(*Lemma de_morgan_4 :  ~(P /\ Q) -> ~P \/ ~Q.
Proof.
 unfold not.
 intros.
 left.
 intro.
 apply H.
 split.
 - assumption.
 - (* loop *) 
*)


Lemma imp_to_not_not_or : (P -> Q) -> ~~(~P \/ Q).
Proof.
  unfold not.
  intros.
  apply H0.
  left.
  intro.
  apply H0.
  right.
  apply H.
  assumption.
Qed.

Print imp_to_not_not_or.

Lemma contraposition : (P -> Q) -> (~Q -> ~P).
Proof.
  unfold not.
  intros.
  apply H0.
  apply H.
  apply H1.
Qed.


(* przykład z wykładu *)
Lemma imp_to_imp_not: (P -> ~Q) -> (~P -> ~Q) -> ~Q.
Proof.
  unfold not.
  intros.
  apply H0.
  - intros.
    apply H; assumption.
  - assumption.
Qed.
 
 (* Now, you may invent and solve your own exercises ! 
    Note that you can trust the tactic tauto: if it fails, then your formula
    is probably not (intuitionnistically) provable *)
(*
Lemma contraposition''' : (~P -> ~Q) <-> (Q -> P).
Proof.
 tauto.
Toplevel input, characters 8-13:
Error: tauto failed.
*)
End propositional_logic.


