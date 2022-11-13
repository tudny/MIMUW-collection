#ifndef SIK_BOMBERMAN_STRING_H
#define SIK_BOMBERMAN_STRING_H


#include "Networkable.h"
#include "Numbers.h"

namespace Models {

    class String : public Networkable {
    public:
        String() = default;

        explicit String(const std::string &value) : value(value) {}

        [[nodiscard]] std::string serialize() const override {
            size_t length = value.length();
            if (length > UINT8_MAX)
                throw SerializationException{};
            UINT8 size{static_cast<uint8_t>(length)};
            return size.serialize() + value;
        }

        void deserialize(Client::Connection &connection) override {
            UINT8 length;
            length.deserialize(connection);
            value = connection.read(length.value);
        }

        std::string value;
    };
}


#endif //SIK_BOMBERMAN_STRING_H
