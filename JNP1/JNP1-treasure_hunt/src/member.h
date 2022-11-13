#ifndef TREASURE_HUNT__MEMBER_H
#define TREASURE_HUNT__MEMBER_H

#include <concepts>
#include <cstdint>

const std::size_t EXPEDITIONS_LIMIT = 25;

template<std::size_t Number>
struct Fibonacci {
    static constexpr uint64_t value = Fibonacci<Number - 1>::value + Fibonacci<Number - 2>::value;
};

template<>
struct Fibonacci<1> {
    static constexpr uint64_t value = 1;
};
template<>
struct Fibonacci<0> {
    static constexpr uint64_t value = 0;
};

template<std::size_t Number, std::size_t Limit>
concept less_then = Number < Limit;

template<treasure_accepted_type ValueType, bool IsArmed>
class Adventurer {
public:
    using strength_t = uint32_t;

    constexpr Adventurer() requires (!IsArmed) {
        _strength = 0;
    }

    constexpr explicit Adventurer(strength_t strength) requires IsArmed {
        _strength = strength;
    }

    constexpr strength_t getStrength() const requires IsArmed {
        return _strength;
    }

    constexpr void loot(SafeTreasure<ValueType> &&treasure) {
        currentLoot += treasure.getLoot();
    }

    constexpr void loot(TrappedTreasure<ValueType> &&treasure) {
        if (_strength != 0) {
            currentLoot += treasure.getLoot();
            _strength >>= 1;
        }
    }

    constexpr ValueType pay() {
        ValueType toPay = currentLoot;
        currentLoot = 0;
        return toPay;
    }

    constexpr static const bool isArmed = IsArmed;

private:
    strength_t _strength;
    ValueType currentLoot = 0;
};

template<treasure_accepted_type ValueType>
using Explorer = Adventurer<ValueType, false>;

template<treasure_accepted_type ValueType, std::size_t CompletedExpeditions>
requires less_then<CompletedExpeditions, EXPEDITIONS_LIMIT>
class Veteran {
public:
    using strength_t = uint32_t;

    constexpr Veteran() = default;

    template<bool IsTrapped>
    constexpr void loot(Treasure<ValueType, IsTrapped> &&treasure) {
        currentLoot += treasure.getLoot();
    }

    constexpr ValueType pay() {
        ValueType toPay = currentLoot;
        currentLoot = 0;
        return toPay;
    }

    constexpr strength_t getStrength() {
        return _strength;
    }

    constexpr static const bool isArmed = true;

private:
    ValueType currentLoot = 0;
    strength_t _strength = Fibonacci<CompletedExpeditions>::value;
};

#endif // TREASURE_HUNT__MEMBER_H