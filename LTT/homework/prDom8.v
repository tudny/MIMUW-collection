(* Aleksander Tudruj - Praca domowa do ćwiczeń 8 *)

Section Classical_tautologies.

Variable DN : forall p, ~~p -> p.


Lemma implication_is_or : forall p q : Prop,
  (~p -> q) -> (p \/ q).
Proof.
unfold not.
intros.
apply DN.
unfold not.
intros.
apply H0.
right.
apply H.
intros.
apply H0.
left.
apply H1.
Qed.

Lemma implication_is_and: forall p q : Prop,
    ~(p -> ~q) -> p /\ q.
Proof.
unfold not.
intros.
apply DN.
unfold not.
intros.
apply H.
intros.
apply H0.
split; assumption.
Qed.

End Classical_tautologies.



Lemma DNEM: (forall p, ~~p -> p) -> (forall p, p \/ ~p).
Proof.
unfold not.
intros.
apply H.
intros.
apply H0.
right.
intros.
apply H0.
left.
apply H1.
Qed.


Lemma EMDN: (forall p, p \/ ~p) -> (forall p, ~~p -> p).
Proof.
unfold not.
intros.
destruct (H p).
- assumption.
- exfalso.
  apply H0.
  apply H1.
Qed.