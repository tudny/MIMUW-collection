#ifndef SIK_BOMBERMAN_PLAYER_H
#define SIK_BOMBERMAN_PLAYER_H

#include <string>
#include <utility>
#include "Position.h"
#include "Networkable.h"
#include "String.h"

namespace Models {

    class Player : public Networkable {
    public:
        Player() = default;

        Player(String name, String address) : name(std::move(name)), address(std::move(address)) {}

        [[nodiscard]] std::string serialize() const override {
            return name.serialize() + address.serialize();
        }

        void deserialize(Client::Connection &connection) override {
            name.deserialize(connection);
            address.deserialize(connection);
        }

        String name, address;
        Position position;
    };

    using PlayerId = UINT8;
    using Score = UINT32;
}


#endif //SIK_BOMBERMAN_PLAYER_H
