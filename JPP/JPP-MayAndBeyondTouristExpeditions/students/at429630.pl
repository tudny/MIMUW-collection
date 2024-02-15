% Aleksander Tudruj

% ==============================================================================
% Solution for task: Zadanie z Prologu (JPP 2023)
%                    Majowe (i nie tylko) wyprawy turystyczne.
% Author:            Aleksander Tudruj (at429630)
% Date:              2023-05-**
% ==============================================================================


% ==============================================================================
%                           Program entrypoint
% ==============================================================================

:- ensure_loaded(library(lists)).

user:runtime_entry(start) :-
    % Read args
    current_prolog_flag(argv, Argv),
    % Validate args
    validate_args(Argv),
    % Load the file
    [TrailsFile|_] = Argv,
    consult(TrailsFile),
    prompt(_Old, ''),
    main_loop,
    halt(0).

% ==============================================================================
%                           Arguments handling
% ==============================================================================

% Program arguments can contain only one argument - the path to a file
validate_args(Argv) :-
    length(Argv, 1),
    !.

validate_args(Argv) :-
    length(Argv, ArgvN), 
    format('Invalid number of arguments: ~w~n', [ArgvN]),
    write('Expected exacly 1 argument.'), nl, nl,
    write('Usage:   ./program <trails-file>'), nl,
    write('Example: ./program tatry.txt'), nl,
    halt(1).

% ==============================================================================
%                           User input handling
% ==============================================================================

% End the program if user enter `koniec.`.
check_if_end(koniec) :-
    !,
    write('Koniec programu. Milych wedrowek!'), nl,
    halt(0).

check_if_end(_).

% Read input line and check for `koniec.` rule.
read_input(Input) :-
    read(Input),
    check_if_end(Input).

% Read trails conditions until user enters correct arguments. 
read_conditions(CondsParsedGood, LenFunGood, XG) :-
    write('Podaj warunki: '),
    read_input(Input),
    (
        % If [user entered parsable arguments]:
        parse_conditons_full(Input, CondsParsed, LenFun, X) -> (
            % Then [unify `read_conditions` arguments]
            CondsParsedGood = CondsParsed,
            LenFunGood = LenFun,
            XG = X
        ) ; (
            % Else [try again]
            read_conditions(CondsParsedGood, LenFunGood, XG)    
        )
    ).

% If given source/destination is `nil` we want to insert ununified variable.
% This will ensure backtrack will take all nodes into consideration.
resolve_nil(X, X) :- X \= nil.
resolve_nil(nil, _).

% Main loop. Reads user input and runs path finding algorithm.
main_loop :-
    write('Podaj miejsce startu: '),
    read_input(Start),
    write('Podaj miejsce koncowe: '),
    read_input(End),
    read_conditions(Types, LenFun, X),
    resolve_nil(Start, ResolvedStart),
    resolve_nil(End, ResolvedEnd),
    do(ResolvedStart, ResolvedEnd, Types, (LenFun, X)),
    main_loop.

% ==============================================================================
%                           Main program logic
% ==============================================================================

do(Start, End, TrailTypesDups, LenFun) :-
    remove_duplicates(TrailTypesDups, TrailTypes), !,
    % Negating `get_all_paths` will ensure all paths are found 
    %  during backtracing, but `do` will return true. 
    \+ get_all_paths(Start, End, TrailTypes, LenFun).


get_all_paths(Start, End, TrailTypes, LenFun) :-
    % We try to find a correct path.
    path_finder(Start, End, TrailTypes, (Length, Path), LenFun),
    LenFun = (LenCondFun, X),
    X = Length,
    (
        % We check is the path length meets the length requirement.
        call(LenCondFun) -> (
            % If yes, print the path.
            nl,
            print_path(Path),
            format('Dlugosc trasy: ~w.~n', [Length]),
            nl
        )
    ),
    % We fail, so Prolog's backtrack will find another path.
    fail.


% ==============================================================================
%                           Conditions
% ==============================================================================

% `dlugosc` must have a proper War and K arguments.
check_len_cond(War, K) :-
    number(K),
    member(War, [eq, lt, gt, le, ge]),
    !.

check_len_cond(War, K) :-
    format(
        'Error: dlugosc ma niepoprawne argumenty - dlugosc(~w,~w)~n', 
        [War, K]
    ),
    !,
    fail.

% `function` checking the trail length. We keep the argument X seperatly,
%   so when we unify it, we can call `trail_length_cond(...)` and check 
%   if it passes.
trail_length_cond(eq, K, X) :- K = X.
trail_length_cond(lt, K, X) :- K > X.
trail_length_cond(gt, K, X) :- K < X.
trail_length_cond(le, K, X) :- K >= X.
trail_length_cond(ge, K, X) :- K =< X.

% This will be injected if path length condition is missing.
always_true(_).

% If no conditions are passed, we use default ones.
parse_conditons_full(nil, [], always_true(X), X).
parse_conditons_full(C, CP, F, X) :-
    parse_conditons(C, CP, F, X).

% A single `rodzaj` condition.
parse_conditons(rodzaj(R), [R], always_true(X), X) :- 
    !.

% A single `dlugosc` condition.
parse_conditons(dlugosc(War, K), [], trail_length_cond(War, K, X), X) :- 
    !, 
    check_len_cond(War, K).


% `rodzaj` condition, but not as the last argument.
parse_conditons((rodzaj(R), RawConditions), [R|Conditions], LenCond, X) :-
    !,
    parse_conditons(RawConditions, Conditions, LenCond, X).

% `dlugosc` conditon is no `dlugosc` was passed after it.
parse_conditons(RConds, Conditions, trail_length_cond(War, K, X), X) :-
    RConds = (dlugosc(War, K), RawConditions),
    check_len_cond(War, K),
    parse_conditons(RawConditions, Conditions, always_true(X), X).

% We inform about multiple `dlugosc` conditons.
parse_conditons((dlugosc(_, _), RawConditions), Conditions, _, X) :-
    parse_conditons(RawConditions, Conditions, trail_length_cond(_, _, X), X),
    format('Error: za duzo warunkow na dlugosc~n', []),
    !,
    fail.

% We check of conditions that are neither `dlugosc` nor `rodzaj`.
parse_conditons((Cond, _), _, _, _) :- wrong_cond(Cond).
parse_conditons(Cond, _, _, _) :- wrong_cond(Cond).

% We need to check if given condition is either `dlugosc` or `rodzaj`
wrong_cond(Cond) :-
    Cond \= dlugosc(_, _),
    Cond \= rodzaj(_),
    !,
    format('Error: niepoprawny warunek - ~w~n', [Cond]),
    !,
    fail.

% Remove duplicates from list. This will make a program **a little** faster.
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


% ==============================================================================
%                            Graph
% ==============================================================================

% Coding `trasa` terms into `edge` terms.
% Coding Ids into pairs will distinguish one-way and two-ways edges.
edge((ID, 0), FROM, TO, WEIGTH, TYPE) :- 
    trasa(ID, FROM, TO, TYPE, jeden, WEIGTH).
edge((ID, 1), FROM, TO, WEIGTH, TYPE) :- 
    trasa(ID, FROM, TO, TYPE, oba, WEIGTH).
edge((ID, 2), FROM, TO, WEIGTH, TYPE) :- 
    trasa(ID, TO, FROM, TYPE, oba, WEIGTH).


% Check if types path meets the requirements.
% If no conditions for type were passes, all edges are ok.
type_meets_conditions(_, []) :- !.

type_meets_conditions(Type, Conditions) :-
    memberchk(Type, Conditions), !.

% Check if list contains at least two elements.
% Path is stored as [N1, E1, N2, E2, ..., Nn],
%   so an empty path is represented as [N1] - a single node without edges.
not_singleton([_, _ | _]).

% Path cannot be empty (check above for definition).
path_finder(Start, End, Conditions, (Length, Path), LenFun) :-
    not_singleton(Path),
    path_with_types(Start, End, Conditions, (Length, Path), LenFun).

% Loading accumulator varaibles for the search.
path_with_types(Start, End, Conds, PathwLen, LenFun) :-
    path_with_types(Start, End, Conds, 0, [], [Start], PathwLen, LenFun).

% End of path. Reversig the result.
path_with_types(End, End, _, Length, _, RPath, (Length, Path), _) :-
    reverse(RPath, Path).

% The `meat` of the algorithm.
path_with_types(Start, End, Conds, LenCur, Visi, RPath, PathwLen, LenFun) :-
    % We check all edges Start->Next
    edge(Id, Start, Next, Weigth, Type),
    % Edge has to have proper type.
    type_meets_conditions(Type, Conds),
    % Edge could not be visited before.
    \+ memberchk(Id, Visi),
    % Updating the length of currently analyzed path.
    Length is LenCur + Weigth,
    % [opt] We can exclude paths, 
    %  that already does not meet the length requirement. 
    limit_too_long_paths(Length, LenFun),
    % Decode the id
    (EdgeId, _) = Id,
    % Run for Next node.
    path_with_types(
        Next, End,
        Conds,
        Length,
        [Id | Visi],
        [Next, (EdgeId, Type) | RPath],
        PathwLen,
        LenFun
    ).

% This should exclude paths that are already too long.
limit_too_long_paths(_, (always_true(X), X)).
limit_too_long_paths(Len, (trail_length_cond(eq, K, X), X)) :- Len =< K.
limit_too_long_paths(Len, (trail_length_cond(lt, K, X), X)) :- Len < K.
limit_too_long_paths(Len, (trail_length_cond(le, K, X), X)) :- Len =< K.
limit_too_long_paths(_, (trail_length_cond(gt, _, X), X)).
limit_too_long_paths(_, (trail_length_cond(ge, _, X), X)).

% ==============================================================================
%                              Print path
% ==============================================================================

% An empty list doesn't even represent an empty path.
print_path([]) :-
    write('This should never happen.'),
    halt(1).

% Single node (at the end of every path)
print_path([A]) :-
    write(A), nl.

% Node with edge -  N->
print_path([N, (Id, T) | R]) :-
    format('~w -(~w,~w)-> ', [N, Id, T]),
    print_path(R).

