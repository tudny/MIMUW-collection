Section praca_domowa.

(* Aleksander Tudruj *)

Variables P Q R A B D X: Prop.


Lemma z1: ~(~P \/ Q) -> ~(P -> Q).
  unfold not.
  intros.
  apply H.
  left.
  intros.
  apply H.
  right.
  apply H0.
  apply H1.
Qed.

Check z1.
Print z1.

Lemma z2: (P -> ~P) -> (~P -> P) -> False.
  unfold not.
  intros.
  apply H.
  - apply H0.
    intros.
    apply H.
    apply H1.
    apply H1.
  - apply H0.
    apply H.
    apply H0.
    intros.
    apply H.
    apply H1.
    apply H1.
Qed.

Check z2.
Print z2.

Lemma z3: ((R -> A) -> P) \/ ((D -> B) -> Q) -> (P -> X) -> ~(A -> X) -> B -> Q.
  unfold not.
  intros.
  destruct H.
  - exfalso.
    apply H1.
    intros.
    apply H0.
    apply H.
    intros. 
    apply H3.
  - apply H.
    intros.
    apply H2.
Qed.

Check z3.
Print z3.

Lemma z4: ~(P -> Q) -> ((Q -> P) -> P ) -> P.
  unfold not.
  intros.
  apply H0.
  intros.
  exfalso.
  apply H.
  intros.
  apply H1.
Qed.

Check z4.
Print z4.

Lemma z5: ~~(~~P-> P).
  unfold not.
  intros H.
  apply H.
  intro.
  exfalso.
  apply H0.
  intros.
  apply H.
  intros.
  apply H1.
Qed.

End praca_domowa.

