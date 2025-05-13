Section pd6.

Variables P Q : Prop.

(* 1. *)

Lemma first1 : P \/ Q -> (~P -> Q).
Proof.
  unfold not.
  intros.
  destruct H.
  - exfalso.
    apply H0.
    apply H.
  - apply H. 
Qed.

Lemma first2 : (~P -> Q) -> P \/ Q.
Proof.
Admitted.

(* 2. *)

Lemma second1 : P /\ Q -> ~(P -> ~Q).
Proof.
  unfold not.
  intros.
  destruct H.
  apply H0.
  - apply H.
  - apply H1. 
Qed.

Lemma second2 : ~(P -> ~Q) -> P /\ Q.
Proof.
Admitted.

(* 3. *)

Lemma third1 : P /\ Q -> (P <-> (P -> Q)).
Proof.
  intros.
  destruct H.
  split.
  - intros.
    apply H0.
  - intros.
    apply H.
Qed.

Lemma third2 : (P <-> (P -> Q)) -> P /\ Q.
Proof.
  intro.
  destruct H.
  split.
  - apply H0.
    intro.
    apply H; apply H1.
  - apply H; apply H0; intro; apply H; apply H1.
Qed.

(* 4. *)

Lemma fouth1 : (P -> Q) -> ((P \/ Q) <-> Q).
Proof.
  intro.
  split.
  - intro.
    destruct H0.
    * apply H.
      apply H0.
    * apply H0.
  - intro.
    right.
    apply H0.
Qed.

Lemma fouth2 : ((P \/ Q) <-> Q) -> (P -> Q).
Proof.
  intros.
  destruct H.
  apply H.
  left.
  apply H0.
Qed.

End pd6.


