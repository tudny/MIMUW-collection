#ifndef SIK_BOMBERMAN_NETWORKABLE_H
#define SIK_BOMBERMAN_NETWORKABLE_H


#include <string>
#include <concepts>
#include "../../client/Connection.h"

namespace Models {

class SerializationException : public std::exception {};
class DeserializationException : public std::exception {};

    class Networkable {
    public:
        [[nodiscard]] virtual std::string serialize() const = 0;
        virtual void deserialize(Client::Connection &connection) = 0;
    };

    template<class T>
    concept networkable = std::derived_from<T, Models::Networkable>;
}


#endif //SIK_BOMBERMAN_NETWORKABLE_H
