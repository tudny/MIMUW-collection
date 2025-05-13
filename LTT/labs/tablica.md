

1. 
```hs
2 K

2 :: (p -> p) -> (p -> p)
2 = \f \x . f (f x)

K :: (p -> q -> p)
K = \x \y . x

2 K -- SIĘ NIE TYPUJE
```

2.
```hs
A = \x . K (K x)

A :: p -> q -> s -> p
```


# Krata (zupełna)

Krata - porządek częściowy `< A, <= >`
sup i inf każdego zbioru dwuelementowego istnieje

Krata jest dystrybuatywna jeśli zachodzi
```
-- wiadomo co
```


## Def. Algebra Boole'a
to systrybuatywna krata B z el. największym i najmniejszym, taka że każdy element `a` ma dopełnienie (oznaczone `-a`).

### Przykłady
```
< P(X), \cup, \cap, \emptyset, X >
< {\emptyset, X}, \cup, \cap, -, \emptyset, X >
< {A \subseteq \N : A skończone lub X-A skończone}, \cup, \cap, -, \emptyset,\N >
```

## Wartościowanie w algebrze Boole'a
`B = <B, \cup, \cap, -, 0, 1>`
```
V : Var -> B
```

Wartości formuły zdaniowej
-- po bożemu


## Praca domowa

Niech $B = < P(X), \cup, \cap, -, \emptyset, X >$
Zad. Dla każdej formuły klasycznego rachunku zdań p
$$
B |= \phi <=> \phi \;\;\text{jest tautologią (w semantyce 0-1)}
$$
