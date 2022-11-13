#include "treasure_hunt.h"
#include <cstdint>

namespace {

constexpr int soloHunt() {
    Explorer<int> a;
    Treasure<int, false> t1(5);
    Treasure<int, false> t2(6);
    SafeTreasure<int> t3(7);
    TrappedTreasure<int> t4(10);

    Encounter<decltype(a), decltype(t1)> e1 = {a, t1};
    Encounter<decltype(t2), decltype(a)> e2 = {t2, a};
    Encounter<decltype(t3), decltype(a)> e3 = {t3, a};
    Encounter<decltype(a), decltype(t4)> e4 = {a, t4};

    expedition(e1, e2, e3, e4);

    return a.pay(); // 18
}

constexpr int holyGrail() {
    TrappedTreasure<int> grail(10000);
    Veteran<int, 6> arthur;
    Adventurer<int, true> morgana(15);

    Encounter<decltype(arthur), decltype(grail)> e1 = {arthur, grail};
    Encounter<decltype(morgana), decltype(grail)> e2 = {morgana, grail};
    Encounter<decltype(arthur), decltype(morgana)> e3 = {arthur, morgana};

    expedition(e2, e1, e1, e3, e3);

    return arthur.pay(); // 10000
}

constexpr int tiringTraps() {
    Adventurer<int, true> a(2);
    TrappedTreasure<int> t1(1);
    TrappedTreasure<int> t2(1);
    TrappedTreasure<int> t3(1);
    TrappedTreasure<int> t4(1);

    Encounter<decltype(a), decltype(t1)> e1 = {a, t1};
    Encounter<decltype(a), decltype(t1)> e2 = {a, t2};
    Encounter<decltype(a), decltype(t1)> e3 = {a, t3};
    Encounter<decltype(a), decltype(t1)> e4 = {a, t4};

    expedition(e1, e2, e3, e4);

    return a.pay(); // 2
}

} // anonimowa przestrzeń nazw

int main() {
    // Poprawne instancje szablonów
    constexpr Treasure<char, false> smallTreasure{'a'};
    constexpr SafeTreasure<int> freeMoney{100};
    constexpr TrappedTreasure<int64_t> bigChest{8000000000};
    constexpr Veteran<int, 18> v;

    // Niepoprawne instancje szablonów
    // constexpr Treasure<float, true> floatingChest{1.0};
    // constexpr Veteran<int, 27> lyingVet{123};

    static_assert(soloHunt() == 18);
    static_assert(holyGrail() == 10000);
    static_assert(tiringTraps() == 2);
}
