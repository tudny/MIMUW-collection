#ifndef SIK_BOMBERMAN_BOMB_H
#define SIK_BOMBERMAN_BOMB_H

#include <utility>

#include "Networkable.h"
#include "Position.h"
#include "Numbers.h"

namespace Models {

    class Bomb : public Networkable {
    public:
        Bomb() = default;

        Bomb(Position position, UINT16 timer) : position(std::move(position)), timer(std::move(timer)) {}

        [[nodiscard]] std::string serialize() const override {
            return position.serialize() + timer.serialize();
        }

        void deserialize(Client::Connection &connection) override {
            position.deserialize(connection);
            timer.deserialize(connection);
        }

        Position position;
        UINT16 timer;
    };

    using BombId = UINT32;
}


#endif //SIK_BOMBERMAN_BOMB_H
