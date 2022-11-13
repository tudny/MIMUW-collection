#ifndef ROBOTS_SERIALIZER_H
#define ROBOTS_SERIALIZER_H

#define PTR_CAST(X, C) ((C *) X)

/* BOOST */
#include <boost/endian/conversion.hpp>

/* STD */
#include <memory>
#include <utility>
#include <vector>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <numeric>

namespace endian = boost::endian;

std::string generate_string_of_uint8(const uint8_t &number);
std::string generate_string_of_uint16(const uint16_t &number);
std::string generate_string_of_uint32(const uint32_t &number);

class DeserializationException : public std::exception {};
class StringLengthException : public std::exception {
public:
    explicit StringLengthException(std::string message) : message(std::move(message)) {}

    ~StringLengthException() noexcept override = default;

    [[nodiscard]] const char *what() const noexcept override {
        return message.c_str();
    }

private:
    std::string message;
};

class Serializable {
public:
    virtual ~Serializable() = default;
    [[nodiscard]] virtual std::string serialized() const = 0;
};

class SerializableString : public Serializable {
public:
    SerializableString() = default;
    explicit SerializableString(std::string str) : _str(std::move(str)) {}

    ~SerializableString() override = default;

    [[nodiscard]] std::string serialized() const override;

    bool operator<(const SerializableString &other) const {
        return _str < other._str;
    }

private:
    std::string _str;
};

class SerializableUINT8 : public Serializable {
public:
    SerializableUINT8() = default;
    explicit SerializableUINT8(const uint8_t &number) : _number(number) {}

    ~SerializableUINT8() override = default;

    [[nodiscard]] std::string serialized() const override;

    bool operator<(const SerializableUINT8 &other) const {
        return _number < other._number;
    }

    uint8_t _number{};
};

class SerializableUINT16 : public Serializable {
public:
    SerializableUINT16() = default;
    explicit SerializableUINT16(const uint16_t &number) : _number(number) {}

    ~SerializableUINT16() override = default;

    [[nodiscard]] std::string serialized() const override;

    bool operator<(const SerializableUINT16 &other) const {
        return _number < other._number;
    }

    uint16_t _number{};
};

class SerializableUINT32 : public Serializable {
public:
    SerializableUINT32() = default;
    explicit SerializableUINT32(const uint32_t &number) : _number(number) {}

    ~SerializableUINT32() override = default;

    [[nodiscard]] std::string serialized() const override;

    bool operator<(const SerializableUINT32 &other) const {
        return _number < other._number;
    }

    SerializableUINT32 operator++() {
        return SerializableUINT32{_number++};
    }

    uint32_t _number{};
};

using SerializablePlayerId = SerializableUINT8;

class SerializablePlayer : public Serializable {
public:
    SerializablePlayer() = default;
    SerializablePlayer(const std::string &name, const std::string &address) :
        _name(name), _address(address) {}
    SerializablePlayer(const SerializableString &name, const SerializableString &address);

    ~SerializablePlayer() override = default;

    [[nodiscard]] std::string serialized() const override;

    SerializableString _name;
    SerializableString _address;
};

class SerializableDirection : public Serializable {
    static const uint8_t DIR_COUNT = 4;

public:
    SerializableDirection() = default;
    explicit SerializableDirection(const uint8_t &dir) : _dir(dir) {
        if (dir >= DIR_COUNT)
            throw DeserializationException{};
    }

    ~SerializableDirection() override = default;

    [[nodiscard]] std::string serialized() const override;

    static SerializableDirection deserialized(const std::string &message);

private:
    uint8_t _dir{};
};

template<class T, class U>
concept Derived = std::is_base_of<U, T>::value;

template<Derived<Serializable> T>
class SerializableList : public Serializable {
public:
    explicit SerializableList(std::vector<T> _list) : list(_list) {}
    SerializableList() = default;
    ~SerializableList() override = default;

    [[nodiscard]] std::string serialized() const override {
        std::string init;
        std::string elements_serialized =
                std::accumulate(list.begin(), list.end(), std::move(init),
                                [](const std::string &acc, const T &element) {
                    return acc + element.serialized();
                });

        return generate_string_of_uint32(static_cast<uint32_t>(list.size())) + elements_serialized;
    }

    std::vector<T> list;
};

template<Derived<Serializable> T, Derived<Serializable> U>
class SerializableMap : public Serializable {
public:
    SerializableMap() = default;
    ~SerializableMap() override = default;

    [[nodiscard]] std::string serialized() const override {
        std::string init;
        std::string elements_serialized =
                std::accumulate(map.begin(), map.end(), std::move(init),
                                [](const std::string& acc, const std::pair<const T, U> &element) {
            return acc + element.first.serialized() + element.second.serialized();
        });

        return generate_string_of_uint32(static_cast<uint32_t>(map.size())) + elements_serialized;
    }

public:
    std::map<T, U> map;
};

using SerializableBombId = SerializableUINT32;

class SerializablePosition : public Serializable {
public:
    SerializablePosition() = default;
    SerializablePosition(const uint16_t &x, const uint16_t &y) : _x(x), _y(y) {}
    SerializablePosition(const SerializableUINT16 &x, const SerializableUINT16 &y);

    ~SerializablePosition() override = default;

    [[nodiscard]] std::string serialized() const override;

    bool operator<(const SerializablePosition &other) const {
        return std::make_pair(_x, _y) < std::make_pair(other._x, other._y);
    }

    SerializableUINT16 _x, _y;
};

class SerializableBomb : public Serializable {
public:
    SerializableBomb() = default;
    SerializableBomb(const SerializablePosition &position, const SerializableUINT16 &timer)
        : _position(position), _timer(timer) {}

    ~SerializableBomb() override = default;

    [[nodiscard]] std::string serialized() const override;

    SerializablePosition _position;
    SerializableUINT16 _timer;
};

using SerializableScore = SerializableUINT32;

/* CLIENT -> GUI */

class DrawMessage : public Serializable {};

class Lobby : public DrawMessage {
public:
    Lobby() = default;

    Lobby(const SerializableString &serverName,
          const SerializableUINT8 &playersCount,
          const SerializableUINT16 &x,
          const SerializableUINT16 &y,
          const SerializableUINT16 &gameLength,
          const SerializableUINT16 &explosionRadius,
          const SerializableUINT16 &bombTimer,
          const SerializableMap<SerializablePlayerId, SerializablePlayer> &players);

    ~Lobby() override = default;

    [[nodiscard]] std::string serialized() const override;

private:
    SerializableString server_name;
    SerializableUINT8 players_count;
    SerializableUINT16 x;
    SerializableUINT16 y;
    SerializableUINT16 game_length;
    SerializableUINT16 explosion_radius;
    SerializableUINT16 bomb_timer;
    SerializableMap<SerializablePlayerId, SerializablePlayer> players;
};

class Game : public DrawMessage {
public:
    Game() = default;

    Game(const SerializableString &serverName,
         const SerializableUINT16 &x,
         const SerializableUINT16 &y,
         const SerializableUINT16 &gameLength,
         const SerializableUINT16 &turn,
         const SerializableMap<SerializablePlayerId, SerializablePlayer> &players,
         const SerializableMap<SerializablePlayerId, SerializablePosition> &playerPositions,
         const SerializableList<SerializablePosition> &blocks,
         const SerializableList<SerializableBomb> &bombs,
         const SerializableList<SerializablePosition> &explosions,
         const SerializableMap<SerializablePlayerId, SerializableScore> &scores);

    ~Game() override = default;

    [[nodiscard]] std::string serialized() const override;

private:
    SerializableString server_name;
    SerializableUINT16 x;
    SerializableUINT16 y;
    SerializableUINT16 game_length;
    SerializableUINT16 turn;
    SerializableMap<SerializablePlayerId, SerializablePlayer> players;
    SerializableMap<SerializablePlayerId, SerializablePosition> player_positions;
    SerializableList<SerializablePosition> blocks;
    SerializableList<SerializableBomb> bombs;
    SerializableList<SerializablePosition> explosions;
    SerializableMap<SerializablePlayerId, SerializableScore> scores;
};

/* KLIENT -> SERVER */

class ClientMessage : public Serializable {};

class JoinClient : public ClientMessage {
    const uint8_t JOIN_MESSAGE_ID = 0;

public:
    JoinClient() = default;
    explicit JoinClient(const std::string &str) : name(str) {}
    ~JoinClient() override = default;

    [[nodiscard]] std::string serialized() const override;

private:
    SerializableString name;
};

class PlaceBombClient : public ClientMessage {
    const uint8_t PLACE_BOMB_MESSAGE_ID = 1;

public:
    PlaceBombClient() = default;
    ~PlaceBombClient() override = default;

    [[nodiscard]] std::string serialized() const override;
};

class PlaceBlockClient : public ClientMessage {
    const uint8_t PLACE_BLOCK_MESSAGE_ID = 2;

public:
    PlaceBlockClient() = default;
    ~PlaceBlockClient() override = default;

    [[nodiscard]] std::string serialized() const override;
};

class MoveClient : public ClientMessage {
    const uint8_t MOVE_MESSAGE_ID = 3;

public:
    MoveClient() = default;
    explicit MoveClient(const SerializableDirection &direction) : direction(direction) {}
    ~MoveClient() override = default;

    [[nodiscard]] std::string serialized() const override;

private:
    SerializableDirection direction;
};

/* GUI -> CLIENT */

struct InputMessage {
    virtual std::shared_ptr<ClientMessage> generate_client_message() = 0;
};

class PlaceBomb : public InputMessage {
public:
    std::shared_ptr<ClientMessage> generate_client_message() override;
};

class PlaceBlock : public InputMessage {
public:
    std::shared_ptr<ClientMessage> generate_client_message() override;
};

class Move : public InputMessage {
public:
    Move() = default;
    explicit Move(const SerializableDirection &_direction) :
        direction(_direction) {}

    std::shared_ptr<ClientMessage> generate_client_message() override;

    SerializableDirection direction;
};

static void assure(bool pred) {
    if (!pred)
        throw DeserializationException{};
}

class InputMessageFactory {

    using input_message_maker_t = const std::function<std::shared_ptr<InputMessage>(const std::string &)>;

public:
    [[nodiscard]] static std::shared_ptr<InputMessage> getInputMessage(const std::string &message) {
        assure(!message.empty());

        static const input_message_maker_t input_message_makers[] = {
                [](const std::string &body) {
                    assure(body.empty());
                    return std::make_shared<PlaceBomb>();
                },
                [](const std::string &body) {
                    assure(body.empty());
                    return std::make_shared<PlaceBlock>();
                },
                [](const std::string &body) {
                    return std::make_shared<Move>(SerializableDirection::deserialized(body));
                }
        };
        static const uint8_t input_message_makers_size =
                sizeof(input_message_makers) / sizeof(input_message_makers[0]);

        auto id = static_cast<uint8_t>(message[0]);
        auto body = message.substr(sizeof(uint8_t));
        assure(id < input_message_makers_size);
        return input_message_makers[id](body);
    }
};

#endif //ROBOTS_SERIALIZER_H
