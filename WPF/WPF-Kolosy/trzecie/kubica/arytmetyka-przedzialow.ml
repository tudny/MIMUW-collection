(* Sygnatura do zadania na laboratorium o arytmetyce przedzia³ów. *)
module type ARYTMETYKA_PRZEDZIALOW = sig
  type wartosc
  val wartosc_dokladnosc: float -> float -> wartosc
  val wartosc_od_do: float -> float -> wartosc
  val wartosc_dokladna: float -> wartosc

  val in_wartosc: wartosc -> float -> bool
  val min_wartosci: wartosc -> float
  val max_wartosci: wartosc -> float
  val sr_wartosci:  wartosc -> float

  val plus: wartosc -> wartosc -> wartosc
  val minus: wartosc -> wartosc -> wartosc
  val razy: wartosc -> wartosc -> wartosc
  val podzielic: wartosc -> wartosc -> wartosc
end;;
