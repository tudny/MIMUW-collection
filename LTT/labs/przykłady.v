Section propositional_logic_examples.

Variables P Q R : Prop.


Lemma id_P: P -> P.
Proof.
intro.
assumption.
Qed.

(* intro odpowiada regule wprowadzenia implikacji *)
(* assumption odpowiada aksjomatowi *)



Check id_P.
Print id_P.

(* fun X:A => X is Coq syntax for lambda X:A.X *)


Lemma imp_elim: (P -> Q) -> P -> Q.
Proof.
intros.
apply H.
assumption.
Qed.


(* intros to wiele razy uzyte intro *)
(* apply H odpowiada regule eliminacji implikacji,
   gdzie H to nazwa dowodu implikacji, którą eliminujemy *) 
(* apply działa dla "dowolnie długich" ciągów implikacji - patrz poniżej *)



Lemma big_imp: (P-> Q-> R) -> P -> Q -> R.
Proof.
intros.
apply H.
- apply H0.
- apply H1.
Qed. 


Lemma or_comm: P \/ Q -> Q \/ P.
Proof.
intros.
destruct H.
 - right. 
   assumption.
 - left.
   assumption.
Qed.
 
(* elim H, gdzie H to nazwa dowodu alternatywy, odpowiada regule 
   eliminacji alternatywy 
   destruct H to elim H; intro.
*)
(* right i left odpowiadają regułom wprowadzenia alternatywy *)


Lemma and_comm: P /\ Q -> Q /\ P.
Proof.
intros.
destruct H.
split.
 - assumption.
 - assumption.
Qed.

(* split odpowiada regule wprowadzenia koniunkcji *)
(* powyżej destruct H, gdzie H to nazwa dowodu koniunkcji, robi to 
   co inna reguła eliminacja koniunkcji z wykładu  *)


(* Ważne: negacja P (w Coqu ~P) to w logice intuicjonistycznej 
   skrót na P -> False ; dlatego trzeba rozwinąć ten skrót taktyką 
   unfold *)

Lemma not_contrad :  ~(P /\ ~P).
Proof.
unfold not.
intro.
destruct H.
apply H0.
assumption.
Qed.


Lemma or_not: (P \/ Q) -> ~P -> Q.
Proof.
intros.
destruct H.
- exfalso.
  apply H0.
  assumption.
- assumption.
Qed.

(* exfalso odpowiada regule eliminacji fałszu, gdy zamiast aktualnego 
   celu dowodowego wolimy udowodnić fałsz *)


(*Można też użyć assert
assert - to wprowadzenie lokalnego lematu, 
         wraz z obowiązkiem udowodnienia go jako kolejnego celu dowodowego
*)

Lemma or_not2: (P \/ Q) -> ~P -> Q.
Proof.
unfold not.
intros.
destruct H.
- assert False.
  * apply H0.
    assumption.
  * destruct H1. 
- assumption.
Qed.

(* destruct H1 powyżej odpowiada eliminacji fałszu, można o tym mysleć
   jako o elim H1 lub o exfalso; assumption *)


(* przykład z wykładu *)

Lemma not_not_exm : ~ ~ (P \/ ~ P).
Proof.
unfold not.
intros.
apply H.
right.
intro.
apply H.
left.
assumption.
Qed.

(* przykład z tablicy z ćwiczeń *)

Lemma dist_and_or_1 : ((P \/ R) /\ (Q \/ R)) -> ((P /\ Q ) \/ R).
Proof.
intro.
destruct H. (* elim. koniunkcji *) 
destruct H. (* elim P \/ R *) 
- destruct H0. (* elim Q \/ R *) 
  * left.
    split.
    assumption.
    assumption.
  * right; assumption.
- right.
  assumption.
Qed.

Lemma dist_and_or_2 : ((P /\ Q ) \/ R)-> ((P \/ R) /\ (Q \/ R)).
Proof.
intros.
split; (destruct H; [destruct H; left; assumption | right; assumption]).
Qed.


End propositional_logic_examples.