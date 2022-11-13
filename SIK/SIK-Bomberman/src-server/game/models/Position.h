#ifndef SIK_BOMBERMAN_POSITION_H
#define SIK_BOMBERMAN_POSITION_H

#include <utility>

#include "Networkable.h"
#include "Numbers.h"

namespace Models {

    class Position : public Networkable {
    public:
        Position() = default;

        Position(UINT16 x, UINT16 y) : x(std::move(x)), y(std::move(y)) {}

        Position(int _x, int _y) {
            assert(x >= 0 && y >= 0);
            x = static_cast<uint16_t>(_x);
            y = static_cast<uint16_t>(_y);
        }

        [[nodiscard]] std::string serialize() const override {
            return x.serialize() + y.serialize();
        }

        void deserialize(Client::Connection &connection) override {
            x.deserialize(connection);
            y.deserialize(connection);
        }

        bool operator<(const Position &other) const {
            auto x_order = x <=> other.x;
            if (x_order == 0) {
                auto y_order = y <=> other.y;
                return y_order < 0;
            }
            return x_order < 0;
        }

        bool operator==(const Position &other) const {
            auto is_less = *this < other;
            auto is_smaller = other < *this;
            return !(is_less || is_smaller);
        }

        UINT16 x, y;
    };
}


#endif //SIK_BOMBERMAN_POSITION_H
