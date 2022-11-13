#ifndef TREASURE_HUNT__TREASURE_H
#define TREASURE_HUNT__TREASURE_H

#include <concepts>

template<typename T>
concept treasure_accepted_type = std::integral<T>;

template<treasure_accepted_type ValueType, bool IsTrapped>
class Treasure {
public:
    constexpr explicit Treasure(ValueType value) {
        _value = value;
    }

    constexpr ValueType evaluate() {
        return _value;
    }

    constexpr ValueType getLoot() {
        ValueType loot = _value;
        _value = 0;
        return loot;
    }

    constexpr static const bool isTrapped = IsTrapped;

private:
    ValueType _value;
};

template<treasure_accepted_type ValueType>
using SafeTreasure = Treasure<ValueType, false>;

template<treasure_accepted_type ValueType>
using TrappedTreasure = Treasure<ValueType, true>;

#endif // TREASURE_HUNT__TREASURE_H