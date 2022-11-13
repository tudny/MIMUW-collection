/* LOCAL */
#include "serializer.h"

#define STRING_OF_INTEGRAL(T, N) char buffer[sizeof(T)]; \
                                 *PTR_CAST(buffer, T) = endian::native_to_big(N); \
                                 return std::string{buffer, sizeof(T)}

std::string generate_string_of_uint8(const uint8_t &number) {
    STRING_OF_INTEGRAL(uint8_t, number);
}

std::string generate_string_of_uint16(const uint16_t &number) {
    STRING_OF_INTEGRAL(uint16_t, number);
}

std::string generate_string_of_uint32(const uint32_t &number) {
    STRING_OF_INTEGRAL(uint32_t, number);
}

std::string SerializableString::serialized() const {
    if (_str.size() > UCHAR_MAX)
        throw StringLengthException{"String is too long. " + std::to_string(_str.size())};

    return generate_string_of_uint8(static_cast<uint8_t>(_str.size())) + _str;
}

std::string SerializableUINT8::serialized() const {
    return generate_string_of_uint8(_number);
}

std::string SerializableUINT16::serialized() const {
    return generate_string_of_uint16(_number);
}

std::string SerializableUINT32::serialized() const {
    return generate_string_of_uint32(_number);
}

SerializablePlayer::SerializablePlayer(
        const SerializableString &name,
        const SerializableString &address) :
        _name(name),
        _address(address) {}

std::string SerializablePlayer::serialized() const {
    auto ser_name = _name.serialized();
    auto ser_address = _address.serialized();

    return ser_name + ser_address;
}

std::string SerializableDirection::serialized() const {
    return generate_string_of_uint8(_dir);
}

SerializableDirection SerializableDirection::deserialized(const std::string &message) {
    if (message.size() != 1)
        throw DeserializationException{};

    return SerializableDirection{static_cast<uint8_t>(message[0])};
}

SerializablePosition::SerializablePosition(const SerializableUINT16 &x, const SerializableUINT16 &y) : _x(x), _y(y) {}

std::string SerializablePosition::serialized() const {
    return _x.serialized() + _y.serialized();
}

std::string SerializableBomb::serialized() const {
    return _position.serialized() + _timer.serialized();
}

Lobby::Lobby(const SerializableString &serverName,
             const SerializableUINT8 &playersCount,
             const SerializableUINT16 &x,
             const SerializableUINT16 &y,
             const SerializableUINT16 &gameLength,
             const SerializableUINT16 &explosionRadius,
             const SerializableUINT16 &bombTimer,
             const SerializableMap<SerializablePlayerId,
             SerializablePlayer> &players) :
                 server_name(serverName),
                 players_count(playersCount),
                 x(x), y(y),
                 game_length(gameLength),
                 explosion_radius(explosionRadius),
                 bomb_timer(bombTimer),
                 players(players) {}

std::string Lobby::serialized() const {
    return generate_string_of_uint8(0) +
           server_name.serialized() +
           players_count.serialized() +
           x.serialized() +
           y.serialized() +
           game_length.serialized() +
           explosion_radius.serialized() +
           bomb_timer.serialized() +
           players.serialized();
}

Game::Game(const SerializableString &serverName,
           const SerializableUINT16 &x,
           const SerializableUINT16 &y,
           const SerializableUINT16 &gameLength,
           const SerializableUINT16 &turn,
           const SerializableMap<SerializablePlayerId,
           SerializablePlayer> &players,
           const SerializableMap<SerializablePlayerId,
           SerializablePosition> &playerPositions,
           const SerializableList<SerializablePosition> &blocks,
           const SerializableList<SerializableBomb> &bombs,
           const SerializableList<SerializablePosition> &explosions,
           const SerializableMap<SerializablePlayerId,
           SerializableScore> &scores) :
                server_name(serverName),
                x(x), y(y),
                game_length(gameLength),
                turn(turn),
                players(players),
                player_positions(playerPositions),
                blocks(blocks),
                bombs(bombs),
                explosions(explosions),
                scores(scores) {}

std::string Game::serialized() const {
    return generate_string_of_uint8(1) +
           server_name.serialized() +
           x.serialized() +
           y.serialized() +
           game_length.serialized() +
           turn.serialized() +
           players.serialized() +
           player_positions.serialized() +
           blocks.serialized() +
           bombs.serialized() +
           explosions.serialized() +
           scores.serialized();
}

std::string JoinClient::serialized() const {
    return generate_string_of_uint8(JOIN_MESSAGE_ID) + name.serialized();
}

std::string PlaceBombClient::serialized() const {
    return generate_string_of_uint8(PLACE_BOMB_MESSAGE_ID);
}

std::string PlaceBlockClient::serialized() const {
    return generate_string_of_uint8(PLACE_BLOCK_MESSAGE_ID);
}

std::string MoveClient::serialized() const {
    return generate_string_of_uint8(MOVE_MESSAGE_ID) + direction.serialized();
}

std::shared_ptr<ClientMessage> PlaceBomb::generate_client_message() {
    return std::make_shared<PlaceBombClient>();
}

std::shared_ptr<ClientMessage> PlaceBlock::generate_client_message() {
    return std::make_shared<PlaceBlockClient>();
}

std::shared_ptr<ClientMessage> Move::generate_client_message() {
    return std::make_shared<MoveClient>(direction);
}
