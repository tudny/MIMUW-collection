#ifndef TREASURE_HUNT__TREASURE_HUNT_H
#define TREASURE_HUNT__TREASURE_HUNT_H

#include "treasure.h"
#include "member.h"
#include <concepts>
#include <utility>

template<typename T>
concept IsTreasure = requires(T x) {
    { Treasure{x} } -> std::same_as<T>;
};

template<typename T>
concept hasStaticIsArmed = requires () {
    { [] () constexpr { return T::isArmed; }() };
};

template<typename T>
concept hasStrength = requires (T x) {
    typename T::strength_t;
    x.getStrength();
};

template<typename T>
concept IsMember = requires(T x) {
    typename T::strength_t;
    requires hasStaticIsArmed<T>;
    { T::isArmed } -> std::convertible_to<bool>;
} && (requires(T x) {
    { x.loot(Treasure<decltype(x.pay()), true>(0)) } -> std::same_as<void>;
} || requires(T x) {
    { x.loot(Treasure<decltype(x.pay()), false>(0)) } -> std::same_as<void>;
});

template<typename T>
concept EncounterSide = IsTreasure<T> || IsMember<T>;

template<typename sideA, typename sideB>
requires EncounterSide<sideA> && EncounterSide<sideB>
class Encounter {
public:
    constexpr Encounter(sideA &a, sideB &b) : a(a), b(b) {}

    sideA &a;
    sideB &b;
};

template<IsTreasure A, IsMember B>
constexpr void run(const Encounter<A, B> &encounter) {
    encounter.b.loot(std::move(encounter.a));
}

template<IsMember A, IsTreasure B>
constexpr void run(const Encounter<A, B> &encounter) {
    encounter.a.loot(std::move(encounter.b));
}

template<IsMember A, IsMember B>
requires hasStrength<A> && hasStrength<B>
constexpr void run(const Encounter<A, B> &encounter) {
    typename A::strength_t strengthA = encounter.a.getStrength();
    typename B::strength_t strengthB = encounter.b.getStrength();

    if (strengthA > strengthB) {
        encounter.a.loot(SafeTreasure<decltype(encounter.b.pay())>(encounter.b.pay()));
    }
    else if (strengthB > strengthA) {
        encounter.b.loot(SafeTreasure<decltype(encounter.a.pay())>(encounter.a.pay()));
    } else {
        return;
    }
}

template<IsMember A, IsMember B>
constexpr void run(const Encounter<A, B> &encounter) {
    if (!encounter.a.isArmed && !encounter.b.isArmed) {
        return;
    }
    else if (encounter.a.isArmed && !encounter.b.isArmed) {
        encounter.a.loot(SafeTreasure<decltype(encounter.b.pay())>(encounter.b.pay()));
    }
    else if (!encounter.a.isArmed && encounter.b.isArmed) {
        encounter.b.loot(SafeTreasure<decltype(encounter.a.pay())>(encounter.a.pay()));
    }
}

template<typename T>
constexpr void expedition(T&& encounter) {
    run(encounter);
}

template<typename T, typename... Args>
constexpr void expedition(T&& encounter, Args& ...args) {
    run(encounter);
    expedition(args...);
}

#endif // TREASURE_HUNT__TREASURE_HUNT_H