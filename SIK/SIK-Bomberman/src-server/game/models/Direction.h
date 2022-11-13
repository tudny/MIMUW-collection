#ifndef SIK_BOMBERMAN_DIRECTION_H
#define SIK_BOMBERMAN_DIRECTION_H


#include "Networkable.h"
#include "Numbers.h"

namespace Models {

    class Direction : public Networkable {
    public:
        [[nodiscard]] std::string serialize() const override {
            assert(false); // should never be used
            return "";
        }

        void deserialize(Client::Connection &connection) override {
            UINT8 _dir;
            _dir.deserialize(connection);
            if (_dir.value < DIR_MIN || _dir.value > DIR_MAX)
                throw DeserializationException{};
            dir = static_cast<Dir>(_dir.value);
        }

        uint8_t DIR_MAX = 3;
        uint8_t DIR_MIN = 0;

        enum Dir {
            Up,
            Right,
            Down,
            Left
        } dir{};
    };
}


#endif //SIK_BOMBERMAN_DIRECTION_H
