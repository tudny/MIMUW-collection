#include "treasure_hunt.h"

#include <cassert>

template<class A, class B>
void assert_left_winner(A a, B b) {
    Treasure<int, false> t1(1);
    Treasure<int, false> t2(2);
    run(Encounter{a, t1});
    run(Encounter{b, t2});
    run(Encounter{a, b});
    assert(a.pay() == 3);
    assert(b.pay() == 0);
    Treasure<int, false> t3(4);
    run(Encounter{b, t3});
    run(Encounter{b, a});
    assert(a.pay() == 4);
    assert(b.pay() == 0);
}

template<class A, class B>
void assert_draw(A a, B b) {
    Treasure<int, false> t1(1);
    Treasure<int, false> t2(2);
    run(Encounter{a, t1});
    run(Encounter{b, t2});
    run(Encounter{a, b});
    assert(a.pay() == 1);
    assert(b.pay() == 2);
}

int main() {
    assert_left_winner<>(Veteran<int, 4>(), Veteran<int, 3>());
    assert_draw<>(Veteran<int, 4>(), Veteran<int, 4>());

    assert_draw<>(Explorer<int>(), Explorer<int>());

    assert_left_winner<>(Adventurer<int, true>(4), Adventurer<int, true>(3));
    assert_draw<>(Adventurer<int, true>(3), Adventurer<int, true>(3));

    assert_left_winner<>(Veteran<int, 0>(), Explorer<int>());
    assert_left_winner<>(Veteran<int, 6>(), Adventurer<int, true>(6));
    assert_left_winner<>(Adventurer<int, true>(100), Veteran<int, 5>());
    assert_draw<>(Veteran<int, 5>(), Adventurer<int, true>(5));

    assert_left_winner<>(Adventurer<int, true>(0), Explorer<int>());
}
