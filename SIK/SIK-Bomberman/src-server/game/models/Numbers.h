#ifndef SIK_BOMBERMAN_NUMBERS_H
#define SIK_BOMBERMAN_NUMBERS_H

#include "Networkable.h"
#include <boost/endian/conversion.hpp>
#include "../../utils/utils.h"

namespace Models {
    template<typename T>
    class UINT : public Networkable {
    public:
        UINT() = default;
        constexpr UINT(const T &value) : value(value) {}

        [[nodiscard]] std::string serialize() const override {
            char buffer[sizeof(T)];
            *((T *) buffer) = boost::endian::native_to_big(value);
            return std::string{buffer, sizeof(T)};
        }

        void deserialize(Client::Connection &connection) override {
            std::string message = connection.read(sizeof(T));
            value = boost::endian::big_to_native(*((T *) message.data()));
        }

        bool operator<(const UINT<T> &other) const {
            return value < other.value;
        }

        auto operator<=>(const UINT<T> &other) const {
            return value <=> other.value;
        }

        UINT<T> operator+() const {
            return *this;
        }

        UINT<T> operator-() const {
            return UINT<T>{-value};
        }

        UINT<T> operator+(const UINT<T> &other) const {
            return UINT<T>{value + other.value};
        }

        UINT<T> operator-(const UINT<T> &other) const {
            return UINT<T>(value - other.value);
        }

        UINT<T> &operator++() {
            ++value;
            return *this;
        }

        UINT<T> &operator--() {
            --value;
            return *this;
        }

        UINT<T> operator++(int) {
            return UINT<T>{value++};
        }

        UINT<T> operator--(int) {
            return UINT<T>{value--};
        }

        UINT<T> &operator+=(const UINT<T> &other) {
            value += other.value;
            return *this;
        }

        UINT<T> &operator-=(const UINT<T> &other) {
            value -= other.value;
            return *this;
        }

        T value;
    };

    using UINT8 = UINT<uint8_t>;
    using UINT16 = UINT<uint16_t>;
    using UINT32 = UINT<uint32_t>;
}


#endif //SIK_BOMBERMAN_NUMBERS_H
