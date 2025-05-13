Section Simple_types.
(** We use a section and a [Variable] so that we do not have to worry
about
    type variables [P Q R S]  *)

Variables P Q R S : Prop. (* rozumiane raczaj jako stałe typowe *)

(* fun X:A => X is Coq syntax for lambda X:A.X *)
(* termy w Coqu są w stylu Churcha - zmienne przy lambdach mają podany
typ *)
Definition K: P -> Q -> P := fun X:P => fun Y:Q => X.

(* typ stałej można napisać, ale Coq umie go sam obliczyć z termu *)
Definition K2 := fun (X:P) (Y:Q) => X.

Goal K = K2.
reflexivity.
Qed.

Definition K3 := fun X:P => fun Y:P => X.

(* Nie uda się udowodnić poniższej równości, bo te dwa termy mają różne
typy
Goal K = K3.
*)


(* Można termy beta redukować *)

Definition KK :=  fun (X: P->P) (Y:Q) => X.

Eval cbv beta in KK (fun x => x).
Eval cbv beta delta in KK (fun x => x).

Goal KK (fun x => x) = fun (_ : Q) (x : P) => x.
reflexivity.


Definition K_1 : P -> Q -> P := K.
Definition K_2 : (Q -> P) -> R -> (Q -> P):= fun X:Q -> P => fun Y:R =>
X.

(* Rekonstrukcja typów - odtwarzanie/"zgadywanie"/przypisanie typów dla
K_2 i K_1 i x.
tak, żeby term t poniżej się typował, czyli ogólnie zadanie:
Dany: term
Znajdź: typ (i środowisko), żeby wszystko do siebie pasowało w danym
lambda rachunku
*)

Definition t := fun x:P => K_2 (K_1 x).
Check t.
Qed.


(* Sprawdzanie typów:
Dane: środowisko, term i typ
Sprawdzenie, że dany term ma dany typ w danym środowisku w danym
rachunku lambda

To właśnie robi Coq gdy podamy taką definicję jak dla K powyżej
Definition K: P -> Q -> P := fun X:P => fun Y:Q => X.
*)


(* Problem inhabitacji to zadanie:
Dany: typ (i środowisko)
Znajdź: term tego typu w danym lambda rachunku
*)


Lemma id_P : P -> P.
Proof.
intro X.
exact X.
(* assumption. *)
Qed.

Print id_P.

(* id_P to nazwa na powyżej zbudowany term typu P-> P *)


Lemma id_PP : (P -> P) -> P -> P.
Proof.
intros.
assumption.
Qed.

Print id_PP.

(* id_PP = fun (_ : P -> P) (H0 : P) => H0
      : (P -> P) -> P -> P *)


Lemma d : (P -> Q) -> P -> Q.
Proof.
intros.
apply H.
assumption.
Qed.

Print d.

Lemma id_PP2 : (P -> P) -> P -> P.
Proof.
intros.
apply H.
assumption.
Qed.

Print id_PP2.

(*id_PP2 = fun (H : P -> P) (H0 : P) => H H0
      : (P -> P) -> P -> P *)

(* to jest inny niż id_PP term typu (P -> P) -> P -> P *)

(* Nie uda się
Goal id_PP = id_PP2.
reflexivity.
*)


Lemma id_PP3 : (P -> P) -> P -> P.
Proof.
intro.
assumption.
Qed.

Print id_PP3.

(*
id_PP3 = fun H : P -> P => H
      : (P -> P) -> P -> P
*)

Goal (fun H : P -> P => H) = (fun (H : P -> P) (H0 : P) => H H0).
reflexivity.
Qed.



(* = zawiera eta równość, ale to nie jest redukcja, którą możemy podać w
Eval *)


Lemma pi1 : P -> P -> P.
Proof.
  intros H J.
  exact H.
Qed.

Print pi1.

(* pi1 = fun H _ : P => H
      : P -> P -> P *)

Lemma pi2 : P -> P -> P.
Proof.
  intros.
  exact H0.
Qed.

Print pi2.

(* pi2 = fun _ H0 : P => H0
      : P -> P -> P *)

(* Nie uda się, to są różne postaci normalne
Goal pi1 = pi2.
reflexivity.
*)



(* ZADANIA NA ĆWICZENIA LUB DOMOWE *)

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