Section IPC2.

Lemma l2: forall P:Prop, 
~~P <-> (exists Q:Prop, P <-> (~Q \/ Q)).
  unfold not.  
  split.
  - intros.
    exists P.
    split.
    * intro.
      right.
      apply H0.
    * intros.
      destruct H0.
      ** exfalso.
         apply H.
         apply H0.
      ** apply H0.
  - intros.
    destruct H.
    destruct H.
    apply H0.
    apply H1.
    left.
    intros.
    apply H0.
    apply H1.
    right.
    apply H2.
Qed.

End IPC2.

Section KodowanieImpredykatywne.

Definition my_False : Prop := forall P:Prop, P.

Definition my_not (P:Prop) := P -> my_False.

Definition my_or (P Q:Prop): Prop := forall R:Prop, 
                                    (P -> R) -> (Q->R) -> R.

Definition my_and (P Q:Prop): Prop := forall R:Prop, 
                                    (P-> Q-> R) -> R.

Definition my_exists (A:Type)(P:A->Prop) : Prop :=
  forall R: Prop, 
    (forall a: A, P a -> R) -> R.

Lemma my_or_intro_l : forall P Q:Prop, P -> my_or P Q.
Proof.
  intros.
  unfold my_or.
  intros.
  apply H0.
  apply H.
Qed.

Lemma my_or_ok : forall P Q:Prop, P \/ Q <-> my_or P Q.
Proof.
  intros.
  split.
  - intros.
    unfold my_or.
    intros.
    destruct H.
    * apply H0.
      apply H.
    * apply H1.
      apply H.
  - unfold my_or.
    intros.
    apply (H (P \/ Q)).
    * intro.
      left.
      apply H0.
    * intro.
      right.
      apply H0.
Qed.

Lemma my_ex_ok :  forall (A:Type)(P:A->Prop),
                   (exists x, P x) <-> (my_exists A P).
Proof.
  intros.
  split.
  - intro.
    destruct H.
    unfold my_exists.
    intros.
    apply (H0 x).
    apply H.
  - unfold my_exists.
    intros.
    apply (H (exists a, P a)).
    intros.
    exists a.
    apply H0.
Qed.

End KodowanieImpredykatywne.



 