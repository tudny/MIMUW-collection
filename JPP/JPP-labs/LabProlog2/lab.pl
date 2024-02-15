
% Napisz predykat enumFromTo(M, N, L) ,który dla ustalonych M,N działa podobnie podobnie jak analogiczna funkcja w Haskellu, n.p.

% ?- enumFromTo(1,0,L).
% L = [].

% ?- enumFromTo(1,8,L).
% L = [1, 2, 3, 4, 5, 6, 7, 8].

enumFromTo(X, Y, []) :- X > Y.
enumFromTo(X, Y, [X|T]) :- X =< Y, Z is X + 1, enumFromTo(Z, Y, T).



% insertionSort(+Lista, ?Posortowana),
% insert(+Lista, +Elem, ?NowaLista)

insert([], X, [X]).
insert([H|T], X, [X|[H|T]]) :- X =< H.
insert([H|T], X, [H|R]) :- X > H, insert(T, X, R).

insertionSort([], []).
insertionSort([H|T], R) :- insertionSort(T, R1), insert(R1, H, R).
