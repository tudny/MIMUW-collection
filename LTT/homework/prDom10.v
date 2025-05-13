(* Aleksander Tudruj - praca domowa do ćwiczeń 10 *)

Section LogikaPierwszegoRzedu.

Variables A: Set.
Variables b: A.
Variables P: A->Prop.
Variables R S: Prop.

(*   Praca domowa   *)

Lemma l1: ~~ (forall a:A, False) -> False.
  unfold not.
  intros.
  apply H.
  intros.
  apply H0.
  apply b.
Qed.

Lemma l2: (~~(forall a, ~~(P a))) -> ~~(P b).
  unfold not.
  intros.
  apply H.
  intros.
  apply (H1 b).
  apply H0.
Qed.

Lemma l3: (~~(forall a, ~~(P a))) -> (forall a, ~~(P a)).
  unfold not.
  intros.
  apply H.
  intros.
  apply (H1 a).
  apply H0.
Qed.

Lemma l4: ~~(exists a:A, False) -> exists a:A, False.
  unfold not.
  intros.
  exfalso.
  apply H.
  intros.
  destruct H0.
  apply H0.
Qed.

Lemma l5: (R -> ~S) -> (~~R -> ~S).
  unfold not.
  intros.
  apply H0.
  intros.
  apply H.
  assumption.
  assumption.
Qed.

Lemma l6: (forall x, (P x -> R)) -> (exists a, P a) -> R.
  intros.
  destruct H0.
  apply (H x).
  apply H0.
Qed.

End LogikaPierwszegoRzedu.
