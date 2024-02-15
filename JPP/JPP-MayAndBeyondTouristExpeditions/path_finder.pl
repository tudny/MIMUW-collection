

trasa(r1, zakopane, brzeziny, rower, oba, 25).
trasa(r2, brzeziny, gasienicowa, rower, oba, 15).
trasa(r3, brzeziny, poroniec, rower, oba, 10).
trasa(r4, poroniec, rusinowa, rower, oba, 6).
trasa(g1, zakopane, kuznice, gorska, oba, 7).
trasa(g2, zakopane, kalatowki, gorska, oba, 5).
trasa(g3, kuznice, gasienicowa, gorska, oba, 7).
trasa(g4, gasienicowa, zawrat, gorska, oba, 6).
trasa(g5, gasienicowa, czarnystaw, gorska, oba, 3).
trasa(g6, zawrat, kozia, gorska, jeden, 5).
trasa(g7, kozia, gasienicowa, gorska, jeden, 7).
trasa(p1, zakopane, gubalowka, piesza, oba, 5).




edge(ID, FROM, TO, WEIGTH, TYPE) :- trasa(ID, FROM, TO, TYPE, jeden, WEIGTH).
edge(ID, FROM, TO, WEIGTH, TYPE) :- trasa(ID, FROM, TO, TYPE, oba, WEIGTH).
edge(ID, FROM, TO, WEIGTH, TYPE) :- trasa(ID, TO, FROM, TYPE, oba, WEIGTH).

edge_finder() :- edge(I, F, T, W, Ty), format('Znaleziono sciezke: ~w ~w ~w ~w ~w ~n', [I, F, T, W, Ty]), fail.


type_meets_conditions(_, []) :- !.

type_meets_conditions(Type, Conditions) :-
    memberchk(Type, Conditions), !.


path_with_types(Start, End, Conditions, Length, Path) :-
    path_with_types(Start, End, Conditions, 0, [Start], Length, Path).

path_with_types(End, End, _, Length, RPath, Length, Path) :-
    reverse(RPath, Path).

path_with_types(Start, End, Conditions, LengthCurr, Visited, FinalLength, Path) :-
    edge(Id, Start, Next, Weigth, Type),
    type_meets_conditions(Type, Conditions),
    \+ memberchk(Next, Visited),
    Length is LengthCurr + Weigth,
    path_with_types(
        Next, End,
        Conditions,
        Length,
        [Next, (Id, Weigth, Type) | Visited],
        FinalLength,
        Path
    ).


bigger_then_10(Length) :- Length > 10.

term_as_argument_test(Condition, Value) :-
    call(Condition, Value).


trail_length_cond(eq, K, X) :- K = X.
trail_length_cond(lt, K, X) :- K < X.
trail_length_cond(gt, K, X) :- K > X.
trail_length_cond(le, K, X) :- K =< X.
trail_length_cond(ge, K, X) :- K >= X.


% Parse input conditions


always_true(_).

parse_conditons(rodzaj(R), [R], always_true(X), X) :- !.

parse_conditons(dlugosc(War, K), [], trail_length_cond(War, K, X), X) :- !.


parse_conditons((rodzaj(R), RawConditions), [R|Conditions], LenCond, X) :-
    parse_conditons(RawConditions, Conditions, LenCond, X), !.

parse_conditons((dlugosc(War, K), RawConditions), Conditions, trail_length_cond(War, K, X), X) :-
    parse_conditons(RawConditions, Conditions, always_true(X), X), !.

parse_conditons((dlugosc(_, _), RawConditions), Conditions, trail_length_cond(_, _, X), X) :-
    parse_conditons(RawConditions, Conditions, _, X),
    format('Error: za duzo warunkow na dlugosc.~n', []),
    !.

parse_conditons((Cond, _), _, _, _) :- wrong_cond(Cond).

parse_conditons(Cond, _, _, _) :- wrong_cond(Cond).


wrong_cond(Cond) :-
    format('Error: niepoprawny warunek - ~w.~n', [Cond]),
    !.

% Remove duplicates from list

remove_duplicates(L1, L2) :-
    remove_duplicates(L1, [], L2).

remove_duplicates([], Acc, Res) :-
    reverse(Acc, Res).

remove_duplicates([H|T], Acc, L2) :-
    memberchk(H, Acc),
    remove_duplicates(T, Acc, L2).

remove_duplicates([H|T], Acc, L2) :-
    \+ memberchk(H, Acc),
    remove_duplicates(T, [H|Acc], L2).


