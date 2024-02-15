

% 1.
% dziecko(Dziecko, Matka, Ojciec)
dziecko(jasio, ewa, jan).
dziecko(stasio, ewa, jan).
dziecko(jan, zofia, adam).
dziecko(adam, anna, jakub).

% Drzewo genealogiczne
%             jasio             stasio
%                  \            /
%                   \          /
%                   ewa   ,  jan
%                             |
%                             |
%                       zofia , adam
%                                 |
%                                 |
%                              anna , jakub


ojciec(Ojciec, Dziecko) :- dziecko(Dziecko, _, Ojciec).

matka(Matka, Dziecko) :- dziecko(Dziecko, Matka, _).

rodzic(Rodzic, Dziecko) :- 
    dziecko(Dziecko, Rodzic, _). 
rodzic(Rodzic, Dziecko) :-
    dziecko(Dziecko, _, Rodzic).

babcia(Dziecko, Babcia) :- 
    matka(Babcia, X), matka(X, Dziecko).
babcia(Dziecko, Babcia) :-
    matka(Babcia, X), ojciec(X, Dziecko).

wnuk(Wnuk, Babcia) :-
    babcia(Wnuk, Babcia).

przodek(Przodek, Potomek) :- rodzic(Przodek, Potomek).
przodek(Przodek, Potomek) :- rodzic(Przodek, X), przodek(X, Potomek).


% 3.

nat(z).
nat(s(X)) :- nat(X).

plus(z, X, X) :- nat(X).
plus(s(X), Y, s(Z)) :- plus(X, Y, Z).

minus(X, Y, Z) :- plus(Y, Z, X).

fib(z, z).
fib(s(z), s(z)).
fib(s(s(K)), N) :- fib(s(K), X), fib(K, Y), plus(X, Y, N).

% 4.

lista([]).
lista([_|T]) :- lista(T).

pierwszy(E, [E|L]) :- lista(L).

ostatni(E, [E]).
ostatni(E, [_|L]) :- ostatni(E, L).

element(E, [E|_]).
element(E, [_|L]) :- element(E, L).

% L1 ++ L2 = L3
scal([], L, L) :- lista(L).
scal([H|L1], L2, L3) :- L3 = [H|L4], scal(L1, L2, L4).

intersect(Z1, Z2) :- element(E, Z1), element(E, Z2).

podziel([], [], []).
podziel([X], [X], []).
podziel([H1|[H2|L]], [H1|NP], [H2|P]) :- podziel(L, NP, P).

prefiks([], L) :- lista(L).
prefiks([H|T], [H|L]) :- prefiks(T, L).

% L = L1 ++ P ++ L2
podlista(P, L) :- scal(T, L2, L), scal(L1, P, T), lista(L1), lista(L2), lista(T).

podciag([], L) :- lista(L).
podciag(P, [_|L]) :- lista(L), podciag(P, L).
podciag([H|P], [H|L]) :- lista(L), lista(P), podciag(P, L).


% wypisz(lista) - wypisuje listę w formacie z przecinkami
wypisz([E]) :- write(E), write('.').
wypisz([H|T]) :- write(H), write(','), wypisz(T).

