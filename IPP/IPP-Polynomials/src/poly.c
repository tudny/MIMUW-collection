/** @file
  Implementacja podstawowych operacji na wielomianach rzadkich wielu zmiennych

  Niezmienniki zachowane w każdym wielomianie to:
  - posortowane jednomiany po potęgach malejąco
  - wielomiany są maksymalnie uproszczone
  - w szczególności nie występują jednomiany o zerowych wielomianach

  @author Aleksander Tudruj <at429620@students.mimuw.edu.pl>
  @date 05/02/2021
 */

#include <stdlib.h>
#include <stdio.h>
#include "poly.h"
#include "memory.h"

/**
 * Sprawdzenie czy wielomian @f$p@f$ ma posortowaną tablicę jednomianów
 * po wykładnikach malejąco. Funkcja przydatna do asercji.
 * @param[in] p : sprawdzany wielomian @f$p@f$
 * @return czy wielomian @f$p@f$ ma dobrze posortowane jednomiany
 * */
static bool isSorted(const Poly *p) {
    if (PolyIsCoeff(p))
        return true;

    for (size_t i = 1; i < p->size; ++i)
        if (MonoGetExp(&p->arr[i - 1]) < MonoGetExp(&p->arr[i]))
            return false;

    return true;
}

/**
 * Sprawdzenie czy jednomian może być skrócony.
 * Funkcja sprawdza czy wprowadzony jednomian może zostać skrócony,
 * to znaczy czy jest postaci @f$p_{i}^0\cdot C@f$, gdzie @f$C\in\mathbb{Z}@f$.
 * @param[in] m : sprawdzany jednomian
 * @return czy jednomian może być skrócony
 * */
static inline bool canMonoBeCut(const Mono *m) {
    return MonoGetExp(m) == 0 && PolyIsCoeff(&m->p);
}

/**
 * Sprawdzenie czy wielomian zawiera dokładnie jeden skracalny jednomian.
 * Funkcja sprawdza czy wielomian zawiera dokładnie jeden jednomian, który
 * spełnia warunek opisany w canMonoBeCut(). Jeżeli tak, można ten wielomian
 * skrócić z postaci @f$P(x_0) = p_{0}^0\cdot C@f$ do @f$P(x_0)=C@f$,
 * gdzie @f$C\in\mathbb{Z}@f$.
 * @param[in] p : sprawdzany wielomian
 * @return czy wielomian może być uproszczony
 * */
static inline bool isPolySingleMonoAndZeroExpCoeff(const Poly *p) {
    if (PolyIsCoeff(p))
        return false;

    return p->size == 1 && canMonoBeCut(&p->arr[0]);
}

/**
 * Sprawdzenie poprawności formy zapisu wielomianu.
 * Sprawdzenie czy przekazany wielomian jest poprawnie zapisany ze względu
 * na zapisane powyżej niezmienniki.
 * Sprawdzenie obejmuje czynności takie jak:
 * - niezerowość tablicy jednomianów
 * - nieskracalność wielomianu
 * - posortowanie jednomianów
 * @param[in] p : sprawdzany wielomian
 * @return czy wielomian jest poprawnie zapisany
 * */
static inline bool hasProperForm(const Poly *p) {
    if (PolyIsCoeff(p))
        return true;

    if (p->size < 1)
        return false;

    if (isPolySingleMonoAndZeroExpCoeff(p))
        return false;

    if (!isSorted(p))
        return false;

    return true;
}

static void printMonoLaTeX(const Mono *m, poly_exp_t idx);
static void printPolyLaTeX(const Poly *p, poly_exp_t idx);

static void printMonoNormalized(const Mono *m, poly_exp_t idx);
static void printPolyNormalized(const Poly *p, poly_exp_t idx);

/**
 * Wypisanie jednomianu LaTeXowo.
 * Wypisanie jednomianu w czytelnej postaci.
 * @param[in] m : wypisywany jednomian.
 * @param[in] idx : identyfikator zmiennej @f$x@f$
 * */
static void printMonoLaTeX(const Mono *m, poly_exp_t idx) {
    printf("x_{%d}^{%d}(", idx, m->exp);
    printPolyLaTeX(&m->p, idx + 1);
    printf(")");
}

/**
 * Wypisanie jednomianu znormalizowanego.
 * Wypisanie jednomianu w czytelnej postaci.
 * @param[in] m : wypisywany jednomian.
 * @param[in] idx : identyfikator zmiennej @f$x@f$
 * */
static void printMonoNormalized(const Mono *m, poly_exp_t idx) {
    printf("(");
    printPolyNormalized(&m->p, idx + 1);
    printf(",%d)", m->exp);
}

/**
 * Wypisanie wielomianu LaTeXowo.
 * Wypisanie wielomianu w czytelnej postaci.
 * @param[in] p : wypisywany wielomian.
 * @param[in] idx : identyfikator zmiennej @f$x@f$
 * */
static void printPolyLaTeX(const Poly *p, poly_exp_t idx) {
    if (PolyIsCoeff(p)) {
        printf("%ld", p->coeff);
    }
    else {
        for (size_t i = 0; i < p->size; ++i) {
            if (i != 0) printf(" + ");
            printMonoLaTeX(&p->arr[i], idx);
        }
    }
}

/**
 * Wypisanie wielomianu znormalizowanego.
 * Wypisanie wielomianu w czytelnej postaci.
 * @param[in] p : wypisywany wielomian.
 * @param[in] idx : identyfikator zmiennej @f$x@f$
 * */
static void printPolyNormalized(const Poly *p, poly_exp_t idx) {
    if (PolyIsCoeff(p)) {
        printf("%ld", p->coeff);
    }
    else {
        // Zmienna musi być przesunięta. W przeciwnym razie wystąpi underflow.
        for (size_t i = p->size; i > 0; --i) {
            printMonoNormalized(&p->arr[i - 1], idx);
            if (i != 1) printf("+");
        }
    }
}

/**
 * Porównanie wykładników.
 * Porównanie wykładników co do wartości. Zwracana wartość @f$x@f$ spełnia
 * - @f$x<0@f$ gdy @f$a<b@f$
 * - @f$x=0@f$ gdy @f$a=b@f$
 * - @f$x>0@f$ gdy @f$a>b@f$
 * @param[in] a : porównywany wykładnik
 * @param[in] b : porównywany wykładnik
 * @return wyżej opisany @f$x@f$
 * */
static int compareExps(poly_exp_t a, poly_exp_t b) {
    if (a < b) return -1;
    return a > b;
}

/**
 * Porównanie jednomianów w porządku malejącym po wykładnikach.
 * Funkcja zwraca wartość @f$x@f$ spełniającą warunek
 * - @f$x<0@f$ gdy @f$a@f$ przed @f$b@f$
 * - @f$x=0@f$ gdy @f$a@f$ na równi z @f$b@f$
 * - @f$x>0@f$ gdy @f$a@f$ za @f$b@f$
 * @param[in] a : wskaźnik na porównywany jednomian
 * @param[in] b : wskaźnik na porównywany jednomian
 * @return wyżej opisany @f$x@f$
 * */
static int compareMonosByExp(const void *a, const void *b) {
    Mono *x = (Mono *)a;
    Mono *y = (Mono *)b;

    return -compareExps(MonoGetExp(x), MonoGetExp(y));
}

/**
 * Sortowanie tablicy jednomianów w porządku zadanym przez compareMonosByExp().
 * @param[in] monos : sortowana tablica
 * @param[in] size : rozmiar sortowanej tablicy
 * */
static void sortMonosByExp(Mono *monos, size_t size) {
    qsort(monos, size, sizeof(Mono), compareMonosByExp);
}

/**
 * Algorytm szybkiego potęgowania.
 * Szybkie wyznaczenie liczby @f$a^n@f$.
 * @param[in] a : podstawa
 * @param[in] n : wykładnik
 * @return @f$a^n@f$
 * */
static inline poly_coeff_t fastPower(poly_coeff_t a, poly_exp_t n) {
    if (n == 0)
        return 1;

    poly_coeff_t score = fastPower(a * a, n / 2);

    if (n % 2 == 1)
        score *= a;

    return score;
}

/**
 * Wyznaczenie maksimum z dwóch liczb.
 * Funkcja zwraca maksimum z dwóch liczb typu poly_exp_t.
 * @param[in] a : pierwsza liczba
 * @param[in] b : druga liczba
 * @return max@f$(a, b)@f$
 * */
static inline poly_exp_t max(poly_exp_t a, poly_exp_t b) {
    return (a > b) ? a : b;
}

static Poly addMonosProperty(size_t count, Mono monos[]);

/**
 * Dodanie dwóch jednomianów stałych.
 * Dla dwóch wielomianów stałych zwracana jest ich suma w postaci wielomianu.
 * Wartości przyjmowane są na własność.
 * @param[in] a : wielomian @f$p(x_0)=c_1@f$
 * @param[in] b : wielomian @f$q(x_0)=c_2@f$
 * @return @f$s(x_0) = c_1 + c_2@f$
 * */
static Poly addPropertyTwoCoeffs(Poly *a, Poly *b) {
    assert(PolyIsCoeff(a) && PolyIsCoeff(b));
    return PolyFromCoeff(a->coeff + b->coeff);
}

/**
 * Dodanie dwóch wielomianów niestałych.
 * Dodanie wielomianu
 * @f$p(x_0) = \sum\limits_{i\in\mathbb{I}} x_0^i\cdot p_{0,i}(x_1)@f$
 * oraz wielomianu
 * @f$q(x_0) = \sum\limits_{i\in\mathbb{I}} x_0^i\cdot q_{0,i}(x_1)@f$.
 * Wartości przyjmowane są na własność.
 * @param[in] a : wielomian @f$p@f$
 * @param[in] b : wielomian @f$q@f$
 * @return @f$p+q@f$
 * */
static Poly addPropertyNonCoeffs(Poly *a, Poly *b) {
    assert(!PolyIsCoeff(a) && !PolyIsCoeff(b));
    size_t monosCnt = a->size + b->size;
    Mono *monos = safeCalloc(monosCnt, sizeof(Mono));

    size_t ptr = 0;
    size_t ptrA = 0;
    size_t ptrB = 0;

    while (ptrA < a->size && ptrB < b->size) {
        if (a->arr[ptrA].exp > b->arr[ptrB].exp)
            monos[ptr++] = a->arr[ptrA++];
        else
            monos[ptr++] = b->arr[ptrB++];
    }

    while (ptrA < a->size)
        monos[ptr++] = a->arr[ptrA++];

    while (ptrB < b->size)
        monos[ptr++] = b->arr[ptrB++];

    Poly res = addMonosProperty(monosCnt, monos);

    safeFree((void **) &monos);
    safeFree((void **) &a->arr);
    safeFree((void **) &b->arr);

    return res;
}

/**
 * Dodanie wielomianu niestałego i stałego.
 * Dodanie wielomianu
 * @f$p(x_0) = \sum\limits_{i\in\mathbb{I}} x_0^i\cdot p_{0,i}(x_1)@f$
 * oraz wielomianu @f$q(x_0)=c@f$, gdzie @f$c\in\mathbb{Z}@f$.
 * Wartości przyjmowane są na własność.
 * @param[in] a : wielomian @f$p@f$
 * @param[in] b : wielomian @f$q@f$
 * @return @f$p+q@f$
 * */
static Poly addPropertyCoeffNonCoeff(Poly *a, Poly *b) {
    assert(PolyIsCoeff(a) && !PolyIsCoeff(b));

    if (PolyIsZero(a))
        return *b;

    Poly tempPoly = { .size = 1, .arr = safeMalloc(sizeof(Mono)) };
    tempPoly.arr[0] = MonoFromPoly(a, 0);

    return addPropertyNonCoeffs(&tempPoly, b);
}

Poly PolyAddProperty(Poly *a, Poly *b) {
    if (PolyIsZero(a))
        return *b;
    else if (PolyIsZero(b))
        return *a;

    if (PolyIsCoeff(a) && PolyIsCoeff(b))
        return addPropertyTwoCoeffs(a, b);
    else if (PolyIsCoeff(a) && !PolyIsCoeff(b))
        return addPropertyCoeffNonCoeff(a, b);
    else if (!PolyIsCoeff(a) && PolyIsCoeff(b))
        return addPropertyCoeffNonCoeff(b, a);
    else
        return addPropertyNonCoeffs(a, b);
}

/**
 * Dodanie jednomiantów w posortowanej tablicy.
 * Dla posortowanej tablicy zwracany jest wielomian będący matematyczną sumą
 * tych jednomianów. Jenomiany NIE mogą być zerowe. Jednomiany w przekazanej
 * tablicy moga ulec zmianie.
 * @param[in] count : liczba jednomianów w tablicy
 * @param[in] monos : tablica jednomianów
 * @return wielomian będący sumą jednomianów
 * */
// already sorted descending, no zero polys
static Poly addMonosProperty(size_t count, Mono monos[]) {
    size_t uniqueExp = count;
    bool *isOk = safeCalloc(count, sizeof(bool));
    for (size_t i = 0; i < count; ++i)
        isOk[i] = true;

    for (size_t i = 1; i < count; ++i) {
        if (isOk[i] && isOk[i - 1]) {
            if (monos[i].exp == monos[i - 1].exp) {
                monos[i].p = PolyAddProperty(&monos[i].p, &monos[i - 1].p);

                isOk[i - 1] = false;
                --uniqueExp;

                if (PolyIsZero(&monos[i].p)) {
                    isOk[i] = false;
                    --uniqueExp;
                }
            }
        }
    }

    Poly res;
    bool isSet = false;

    if (uniqueExp == 0) {
        isSet = true;
        res = PolyZero();
    }

    if (uniqueExp == 1)
        for (size_t i = 0; i < count; ++i)
            if (monos[i].exp == 0 && canMonoBeCut(&monos[i]) && isOk[i])
                isSet = true, res = monos[i].p;

    if (!isSet) {
        res = (Poly) {.size = uniqueExp, .arr = safeCalloc(uniqueExp,
                                                           sizeof(Mono))};
        size_t ptr = 0;

        for (size_t i = 0; i < count; ++i)
            if (isOk[i])
                res.arr[ptr++] = monos[i];
    }

    safeFree((void **) &isOk);

    return res;
}

/**
 * Identyczność jednomianów.
 * @param[in] m - jednomian
 * @return jednomian [m]
 * */
static Mono monoIdentity(Mono m) {
    return m;
}

/**
 * Kopia jednomianu.
 * @param[in] m - jednomian
 * @return kopię jednomianu [m]
 * */
static Mono monoDeepClone(Mono m) {
    return MonoClone(&m);
}

/**
 * Dodanie jednomiantów z tablicy.
 * Dla tablicy zwracany jest wielomian będący matematyczną sumą
 * tych jednomianów. Jenomiany MOGĄ być zerowe. Jednomiany w przekazanej
 * tablicy moga ulec zmianie. Jeżeli tablica nie jest posortowana należy ustawić
 * flagę sort na wartość true, aby funkcja zadziałała poprawnie.
 * @param[in] count : liczba jednomianów w tablicy
 * @param[in] monos : tablica jednomianów
 * @param[in] sort  : czy tablica ma zostać posortowana
 * @return wielomian będący sumą jednomianów
 * */
static Poly polyAddMonosPropertySort(size_t count, Mono *monos, bool sort) {
    size_t countCpy = 0;

    for (size_t i = 0; i < count; ++i)
        if (!PolyIsZero(&monos[i].p))
            monos[countCpy++] = monos[i];

    if (sort)
        sortMonosByExp(monos, countCpy);

    Poly res = addMonosProperty(countCpy, monos);
    return res;
}

/**
 * Dodanie jednomiantów z tablicy.
 * Dla tablicy zwracany jest wielomian będący matematyczną sumą
 * tych jednomianów. Jenomiany MOGĄ być zerowe. Jednomiany w przekazanej
 * tablicy moga ulec zmianie. Jeżeli tablica nie jest posortowana należy ustawić
 * flagę sort na wartość true, aby funkcja zadziałała poprawnie.
 * @param[in] count : liczba jednomianów w tablicy
 * @param[in] monos : tablica jednomianów
 * @param[in] sort  : czy tablica ma zostać posortowana
 * @param[in] f     : funkcja przekształcająca jednomian
 * @return wielomian będący sumą jednomianów
 * */
static Poly polyAddMonosOptSort(size_t count, const Mono monos[], bool sort, Mono (*f)(Mono)) {
    Mono *monosCpy = safeCalloc(count, sizeof(Mono));

    for (size_t i = 0; i < count; ++i)
        monosCpy[i] = (*f)(monos[i]);

    Poly res = polyAddMonosPropertySort(count, monosCpy, sort);

    safeFree((void **) &monosCpy);
    return res;
}

/**
 * Pomnożenie wielomianu przez stałą.
 * Funkcja mnoży wielomian przez stałą zachowując przy tym konstrolę
 * nad przekręceniem się liczby całkowitej i normalizacją wielomianu.
 * Funkcja przyjmuje wielomian na wałasność.
 * @param[in] p : wielomian @f$p@f$
 * @param[in] c : mnożona stała
 * @return @f$p\cdot c@f$
 * */
static Poly multConstProperty(Poly *p, poly_coeff_t c) {
    assert(hasProperForm(p));

    if (c == 0) {
        PolyDestroy(p);
        return PolyZero();
    }

    if (PolyIsCoeff(p))
        return PolyFromCoeff(p->coeff * c);

    for (size_t i = 0; i < p->size; i++)
        p->arr[i].p = multConstProperty(&p->arr[i].p, c);

    Poly res = polyAddMonosOptSort(p->size, p->arr, false, &monoIdentity);
    safeFree((void **) &p->arr);
    return res;
}

/**
 * Mnożenie dwóch wielomianów stałych.
 * @param[in] p : wielomian stały @f$p(x_0)=c_1@f$
 * @param[in] q : wielomian stały @f$q(x_0)=c_2@f$
 * @return @f$p\cdot q@f$
 * */
static Poly mulCoeffPoly(const Poly *p, const Poly *q) {
    assert(PolyIsCoeff(p) && PolyIsCoeff(q));
    return PolyFromCoeff(p->coeff * q->coeff);
}

/**
 * Mnożenie wielomianu stałego i niestałego.
 * @param[in] p : wielomian niestały @f$p@f$
 * @param[in] q : wielomian stały @f$q@f$
 * @return @f$p\cdot q@f$
 * */
static Poly mulNonCoeffAndCoeffPoly(const Poly *p, const Poly *q) {
    assert(!PolyIsCoeff(p) && PolyIsCoeff(q));
    Poly toMultiply = PolyClone(p);
    return multConstProperty(&toMultiply, q->coeff);
}

/**
 * Mnożenie dwóch wielomianów niestałych.
 * @param[in] p : wielomian niestały @f$p@f$
 * @param[in] q : wielomian niestały @f$q@f$
 * @return @f$p\cdot q@f$
 * */
static Poly mulTwoNonCoeffPoly(const Poly *p, const Poly *q) {
    assert(!PolyIsCoeff(p) && !PolyIsCoeff(q));
    assert(isSorted(p) && isSorted(q));

    size_t ws = 0, allSize = p->size * q->size;
    Mono *allMonos = safeCalloc(allSize, sizeof(Mono));

    for (size_t i = 0; i < p->size; i++) {
        for (size_t j = 0; j < q->size; j++) {
            allMonos[ws].p = PolyMul(&p->arr[i].p, &q->arr[j].p);
            allMonos[ws].exp = p->arr[i].exp + q->arr[j].exp;
            ws++;
        }
    }

    Poly s = PolyAddMonos(allSize, allMonos);

    safeFree((void **) &allMonos);
    return s;
}

/**
 * Pomocnicza funkcja do obliczania stopnia wielomianu.
 * Funkcja wywołuje się rekurencyjnie, aż do odpowiedniego poziomu,
 * gdzie poprawia wynik zapisany w zmiennej.
 * @param[in] p : wielomian @f$p@f$, którego stopień sprawdzamy
 * @param[in] actIdx : aktualna zmienna @f$x_{actIdx}@f$
 * @param[in] varIdx : szukana zmienna
 * @param[out] acc : wskaźnik na zmienną z wynikiem
 * @return stopień wielomianu @f$p@f$ po zmiennej @f$x_{varIdx}@f$
 * */
static void degBy(const Poly *p, size_t actIdx, size_t varIdx, poly_exp_t *acc) {
    assert(hasProperForm(p));

    if (PolyIsCoeff(p))
        return;

    if (actIdx == varIdx) {
        // jednomiany są posortowane, więc pierwszy na największy wykładnik
        *acc = max(*acc, MonoGetExp(&p->arr[0]));
    }
    else {
        for (size_t i = 0; i < p->size; ++i) {
            degBy(&p->arr[i].p, actIdx + 1, varIdx, acc);
        }
    }
}

void PolyDestroy(Poly *p) {
    if (PolyIsCoeff(p))
        return;

    for (size_t i = 0; i < p->size; ++i)
        MonoDestroy(&p->arr[i]);

    safeFree((void **) &p->arr);
}

Poly PolyClone(const Poly *p) {
    if (PolyIsCoeff(p))
        return PolyFromCoeff(p->coeff);

    Poly res = {.size = p->size, .arr = safeCalloc(p->size, sizeof(Mono))};
    for (size_t i = 0; i < p->size; ++i) {
        res.arr[i] = MonoClone(&p->arr[i]);
    }

    return res;
}

// Niech p, q będą miały posortowane tablice po współczynnikach malejąco.
Poly PolyAdd(const Poly *p, const Poly *q) {
    Poly a = PolyClone(p);
    Poly b = PolyClone(q);

    return PolyAddProperty(&a, &b);
}

Poly PolyAddMonos(size_t count, const Mono monos[]) {
    return polyAddMonosOptSort(count, monos, true, &monoIdentity);
}

Poly PolyMul(const Poly *p, const Poly *q) {
    assert(hasProperForm(p) && hasProperForm(q));

    if (PolyIsCoeff(p) && PolyIsCoeff(q)) {
        return mulCoeffPoly(p, q);
    }
    else if (PolyIsCoeff(p) && !PolyIsCoeff(q)) {
        return mulNonCoeffAndCoeffPoly(q, p);
    }
    else if (PolyIsCoeff(q) && !PolyIsCoeff(p)) {
        return mulNonCoeffAndCoeffPoly(p, q);
    }
    else {
        return mulTwoNonCoeffPoly(p, q);
    }
}

Poly PolyNeg(const Poly *p) {
    assert(hasProperForm(p));
    Poly a = PolyClone(p);
    return PolyNegProperty(&a);
}

Poly PolyNegProperty(Poly *p) {
    return multConstProperty(p, -1);
}

Poly PolySub(const Poly *p, const Poly *q) {
    Poly qq = PolyNeg(q);
    Poly res = PolyAdd(p, &qq);
    PolyDestroy(&qq);
    return res;
}

Poly PolySubProperty(Poly *p, Poly *q) {
    Poly neg = PolyNegProperty(q);
    return PolyAddProperty(p, &neg);
}

poly_exp_t PolyDegBy(const Poly *p, size_t varIdx) {
    assert(hasProperForm(p));

    if (PolyIsZero(p))
        return -1;

    poly_exp_t cnt = 0;
    degBy(p, 0, varIdx, &cnt);

    return cnt;
}

poly_exp_t PolyDeg(const Poly *p) {
    assert(hasProperForm(p));

    if (PolyIsZero(p))
        return -1;

    if (PolyIsCoeff(p))
        return 0;

    poly_exp_t cnt = 0;
    for (size_t i = 0; i < p->size; ++i) {
        poly_exp_t deg = PolyDeg(&p->arr[i].p) + MonoGetExp(&p->arr[i]);
        cnt = max(cnt, deg);
    }

    return cnt;
}

bool PolyIsEq(const Poly *p, const Poly *q) {
    assert(hasProperForm(p) && hasProperForm(q));

    if (PolyIsCoeff(p) ^ PolyIsCoeff(q))
        return false;

    if (PolyIsCoeff(p) && PolyIsCoeff(q)) {
        return p->coeff == q->coeff;
    }
    else {
        if (p->size != q->size)
            return false;

        for (size_t i = 0; i < p->size; i++) {
            if (MonoGetExp(&p->arr[i]) != MonoGetExp(&q->arr[i]))
                return false;

            if (!PolyIsEq(&p->arr[i].p, &q->arr[i].p))
                return false;
        }
    }

    return true;
}

Poly PolyAt(const Poly *p, poly_coeff_t x) {
    assert(hasProperForm(p));

    if (PolyIsCoeff(p))
        return PolyClone(p);

    Poly res = PolyZero();
    for (size_t i = 0; i < p->size; ++i) {
        poly_exp_t exp = MonoGetExp(&p->arr[i]);
        poly_coeff_t substitution = fastPower(x, exp);

        Poly toSubstitute = PolyFromCoeff(substitution);
        Poly multiplied = PolyMul(&p->arr[i].p, &toSubstitute);
        Poly temp = PolyAdd(&res, &multiplied);

        PolyDestroy(&multiplied);
        PolyDestroy(&res);
        res = temp;
    }

    return res;
}

void PrintPolyLaTeX(const Poly *p, char *label) {
    printf("[%s]: ", label);
    printPolyLaTeX(p, 0);
    printf("\n");
}

void PrintPolyNormalized(const Poly *p) {
    printPolyNormalized(p, 0);
}

Poly PolyCloneMonos(size_t count, const Mono monos[]) {
    if (monos == NULL || count == 0)
        return PolyZero();

    return polyAddMonosOptSort(count, monos, true, &monoDeepClone);
}

Poly PolyOwnMonos(size_t count, Mono *monos) {
    if (monos == NULL)
        return PolyZero();
    if (count == 0) {
        safeFree((void **) &monos);
        return PolyZero();
    }

    Poly res = polyAddMonosPropertySort(count, monos, true);
    safeFree((void **) &monos);
    return res;
}

// COMPOSE module

/** Maksymalna potęga dwójki dla liczb z zakresów. */
static const poly_exp_t MAX_POWER_OF_TWO = 30;

/**
 * Szukanie największej potęgi dwójki, mniejszej równej niż [a].
 * @param[in] a - oraniczenie na potęgę
 * @return potęga dwójki mniejsza równa [a]
 * */
static poly_exp_t lastLeqPowerOfTwo(poly_exp_t a) {
    for (poly_exp_t k = MAX_POWER_OF_TWO; k > 0; k--) {
        if ((1 << k) & a)
            return k;
    }

    return 0;
}

/**
 * Wyliczenie potęg wielomianu do potęg dwójki.
 * Obliczenie potęg wielomianu [p] postaci @f$p^(2^i)@f$. dla @f$i<=n@f$.
 * @param[in] p - potęgowany wielomian
 * @param[in] n - maksymalne @f$i@f$
 * */
static Poly* polyPowersOfTwo(Poly *p, poly_exp_t n) {
    Poly *powers = safeCalloc(n + 1, sizeof(Poly));

    powers[0] = PolyClone(p);

    for (poly_exp_t i = 1; i <= n; ++i) {
        powers[i] = PolyMul(&powers[i - 1], &powers[i - 1]);
    }

    return powers;
}

/**
 * Wyczyszczenie tablicy potęg wielomianu.
 * @param[in] count - liczba potęg
 * @param[in] powers - potęgi
 * */
static void clearPowers(size_t count, Poly *powers) {
    for (size_t i = 0; i <= count; ++i) {
        PolyDestroy(&powers[i]);
    }

    safeFree((void **) &powers);
}

/**
 * Szybkie potęgowanie wielomianu.
 * Szybkie potęgowanie wielomianu znając tablicę potęg danego wielomianu do
 * potęg dwójki.
 * @param[in] polyBiPowers - tablica z potęgi wielomianu do potęg dwójki
 * @param[in] exponent - potęga, do której podnosimy
 * @param[in] maxPowerOfTwo - maksymalna znana potęga dwójki
 * @return wielomian podsniesiony do potęgi [exponent]
 * */
static Poly fastPowerPoly(Poly *polyBiPowers,
                          poly_coeff_t exponent,
                          poly_coeff_t maxPowerOfTwo) {
    Poly result = PolyFromCoeff(1);

    for (poly_coeff_t k = 0; k <= maxPowerOfTwo; ++k) {
        if (exponent & (1 << k)) {
            Poly newRes = PolyMul(&result, &polyBiPowers[k]);
            PolyDestroy(&result);
            result = newRes;
        }
    }

    return result;
}

/**
 * Podstawienie wielomianów [substitutes] pod kolejne zmienne wielomianu [base].
 * Podstawienie kolejnych [depth] wielomianów z tablicy [substitutes] pod kolejne
 * zmienne x_i wielomianu [base]. W przypadku zbyt małej liczby wielomianów
 * w tablicy [substitutes], wstawiane są wielomiany zerowe.
 * @param[in] base : wielomian, do którego podstawiamy
 * @param[in] depth : liczba podstawianych wielomianów
 * @param[in] substitutes : tablica podstawianych wielomianów
 * @return wielomian po podstawieniu
 * */
static Poly polyCompose(const Poly *base, size_t depth, const Poly *substitutes) {
    if (PolyIsCoeff(base))
        return PolyClone(base);

    Poly substitute = (depth == 0) ? PolyZero() : substitutes[0];

    Poly res = PolyZero();

    poly_exp_t maxPower = 0;
    for (size_t k = 0; k < base->size; k++) {
        maxPower = max(base->arr[k].exp, maxPower);
    }

    poly_exp_t maxPowerOfTwo = lastLeqPowerOfTwo(maxPower);
    Poly *substitutePowers = polyPowersOfTwo(&substitute, maxPowerOfTwo);

    for (size_t k = 0; k < base->size; k++) {
        Poly xPower = fastPowerPoly(substitutePowers,
                                    base->arr[k].exp,
                                    maxPowerOfTwo);
        Poly sub = polyCompose(&base->arr[k].p,
                               depth == 0 ? depth : depth - 1,
                               substitutes + 1);
        Poly new = PolyMul(&sub, &xPower);
        PolyDestroy(&xPower);
        PolyDestroy(&sub);

        res = PolyAddProperty(&res, &new);
    }

    clearPowers(maxPowerOfTwo, substitutePowers);

    return res;
}

Poly PolyCompose(const Poly *p, size_t k, const Poly q[]) {
    return polyCompose(p, k, q);
}
