#ifndef SIK_BOMBERMAN_CONTAINERS_H
#define SIK_BOMBERMAN_CONTAINERS_H

#include "Networkable.h"

namespace Models {

    template<networkable K, networkable V>
    class Map : public Networkable {
    public:
        [[nodiscard]] std::string serialize() const override {
            std::string result{UINT32{static_cast<uint32_t>(map.size())}.serialize()};
            for (const auto &[key, value] : map) {
                result += key.serialize() + value.serialize();
            }
            return result;
        }

        void deserialize(Client::Connection &connection) override {
            UINT32 length;
            length.deserialize(connection);
            for (UINT32 i = 0; i < length; ++i) {
                K key; V value;
                key.deserialize(connection);
                value.deserialize(connection);
                map[key] = value;
            }
        }

        std::map<K, V> map;
    };

    template<typename T>
    class PtrList : public Networkable {
    public:
        [[nodiscard]] std::string serialize() const override {
            std::string result{UINT32{static_cast<uint32_t>(list.size())}.serialize()};
            for (const auto &element : list) {
                result += element->serialize();
            }
            return result;
        }

        void deserialize(Client::Connection &connection) override {
            // Should never be used.
            (void) connection;
            assert(false);
        }

        std::vector<std::shared_ptr<T>> list;
    };

    template<networkable T>
    class List : public Networkable {
    public:
        [[nodiscard]] std::string serialize() const override {
            std::string result{UINT32{static_cast<uint32_t>(list.size())}.serialize()};
            for (const auto &element : list) {
                result += element.serialize();
            }
            return result;
        }

        void deserialize(Client::Connection &connection) override {
            UINT32 length;
            length.deserialize(connection);
            for (UINT32 i = 0; i < length; ++i) {
                T element;
                element.deserialize(connection);
            }
        }

        std::vector<T> list;
    };
}


#endif //SIK_BOMBERMAN_CONTAINERS_H
