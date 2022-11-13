/* BOOST */
#include <boost/format.hpp>

/* LOCAL */
#include "GameState.h"
#include "../utils/utils.h"

GameState::GameState() = default;

void Event::handle(GameState &game_state) {
    e_handler(game_state);

    game_state.is_lobby = false;
    game_state.do_send = true;
}

BombPlaced::BombPlaced(const SerializableBombId& id,
                       const SerializablePosition &position) :
                            id(id),
                            position(position) {}

void BombPlaced::e_handler(GameState &game_state) {
    LOG(debug) << "[SERVER->] Received BOMB PLACED.";

    game_state.bombs.map[id] = SerializableBomb{position, game_state.bomb_timer};
}

BombExploded::BombExploded(
        const SerializableBombId& id,
        const SerializableList<SerializablePlayerId> &robotsDestroyed,
        const SerializableList<SerializablePosition> &blocksDestroyed) :
            id(id),
            robots_destroyed(robotsDestroyed),
            blocks_destroyed(blocksDestroyed) {}

/**
 * This function calculates all positions of explosion and adds them to the game state.
 * @param game_state the state of game to be updated
 * @param explosion_position the position of a explosion
 */
static void calculate_explosion(GameState &game_state, SerializablePosition &explosion_position) {
    // We generically iterate through all directions.
    static int dx[] = { 0, -1, 0, 1 };
    static int dy[] = { -1, 0, 1, 0 };
    static int k = sizeof(dx) / sizeof(dx[0]);

    // int is big enough to handle both uint16_t and -1 * uint16_t
    for (int i = 0; i < k; ++i) {
        int x = static_cast<int>(explosion_position._x._number);
        int y = static_cast<int>(explosion_position._y._number);

        for (int j = 0; j <= game_state.explosion_radius._number; ++j) {
            if (x < 0 || y < 0 || x >= game_state.size_x._number || y >= game_state.size_y._number)
                break;

            auto current_serializable = SerializablePosition{
                    static_cast<uint16_t>(x), static_cast<uint16_t>(y)
            };

            game_state.explosions.insert(current_serializable);
            if (game_state.blocks.contains(current_serializable))
                break;

            x += dx[i];
            y += dy[i];
        }
    }
}

void BombExploded::e_handler(GameState &game_state) {
    LOG(debug) << "[SERVER->] Received BOMB EXPLODED.";

    auto exploding_bomb = game_state.bombs.map.find(id);
    game_state.explosions.insert(exploding_bomb->second._position);
    game_state.bombs.map.erase(exploding_bomb);

    for (auto &_id : robots_destroyed.list) {
        if (game_state.players_destroyed_in_turn.insert(_id).second) {
            ++game_state.scores.map[_id];
        }
    }

    for (auto &position : blocks_destroyed.list) {
        game_state.blocks_to_remove.insert(position);
        game_state.explosions.insert(position);
    }

    calculate_explosion(game_state, exploding_bomb->second._position);
}

PlayerMoved::PlayerMoved(
        const SerializablePlayerId& id,
        const SerializablePosition &position) :
            id(id),
            position(position) {}

void PlayerMoved::e_handler(GameState &game_state) {
    LOG(debug) << "[SERVER->] Received PLAYER MOVED.";

    game_state.player_position.map[id] = position;
}

BlockPlaced::BlockPlaced(const SerializablePosition &position) : position(position) {}

void BlockPlaced::e_handler(GameState &game_state) {
    LOG(debug) << "[SERVER->] Received BLOCK PLACED.";

    game_state.blocks.insert(position);
}


/* SERVER MESSAGE */

Hello::Hello(const SerializableString &serverName,
             const SerializableUINT8 &playersCount,
             const SerializableUINT16 &sizeX,
             const SerializableUINT16 &sizeY,
             const SerializableUINT16 &gameLength,
             const SerializableUINT16 &explosionRadius,
             const SerializableUINT16 &bombTimer) :
        server_name(serverName),
        players_count(playersCount),
        size_x(sizeX),
        size_y(sizeY),
        game_length(gameLength),
        explosion_radius(explosionRadius),
        bomb_timer(bombTimer) {}

void Hello::handle(GameState &game_state) {
    LOG(debug) << "[SERVER->] Received HELLO.";

    game_state.server_name = server_name;
    game_state.players_count = players_count;
    game_state.size_x = size_x;
    game_state.size_y = size_y;
    game_state.game_length = game_length;
    game_state.explosion_radius = explosion_radius;
    game_state.bomb_timer = bomb_timer;

    game_state.is_lobby = true;
    game_state.do_send = true;
}

AcceptedPlayer::AcceptedPlayer(
        const SerializablePlayerId& id,
        const SerializablePlayer& player) :
            id(id),
            player(player) {}

void AcceptedPlayer::handle(GameState &game_state) {
    LOG(debug) << "[SERVER->] Received ACCEPTED PLAYER.";

    game_state.players.map[id] = player;
    game_state.scores.map[id]; // set score to 0
    game_state.player_position.map[id] = { 0, 0 };
}

GameStarted::GameStarted(
        const SerializableMap<SerializablePlayerId,
        SerializablePlayer> &players) :
            players(players) {}

void GameStarted::handle(GameState &game_state) {
    LOG(debug) << "[SERVER->] Received GAME STARTED.";

    game_state.players = players;
    game_state.scores.map.clear();
    game_state.player_position.map.clear();
    game_state.bombs.map.clear();

    for (auto &[id, player] : game_state.players.map) {
        game_state.scores.map[id]; // set score to 0
        game_state.player_position.map[id] = { 0, 0 };
    }

    game_state.is_lobby = false;
    game_state.do_send = false;
}

Turn::Turn(
        const SerializableUINT16 &turn,
        const List<EventPtr> &events) :
            turn(turn),
            events(events) {}

void Turn::handle(GameState &game_state) {
    LOG(debug) << "[SERVER->] Received TURN.";

    game_state.turn = turn;

    game_state.explosions.clear();
    game_state.players_destroyed_in_turn.clear();

    for (auto &[bomb_id, bomb] : game_state.bombs.map) {
        --bomb._timer._number;
    }

    for (auto &event : events) {
        event->handle(game_state);
    }

    for (auto &position : game_state.blocks_to_remove) {
        game_state.blocks.erase(position);
    }
    game_state.blocks_to_remove.clear();
    game_state.players_destroyed_in_turn.clear();
    game_state.is_lobby = false;
}

GameEnded::GameEnded(
        const SerializableMap<SerializablePlayerId, SerializableScore> &scores) :
            scores(scores) {}

void GameEnded::handle(GameState &game_state) {
    LOG(debug) << "[SERVER->] Received GAME ENDED.";

    game_state.scores = scores;
    game_state.is_lobby = true;
    game_state.players.map.clear();
}

/* DESERIALIZATION */

#define INTEGRAL_OF_STRING(S, T) const char *buffer = (S).data();                   \
                                 if ((S).size() != sizeof(T))                       \
                                    throw DeserializationException{};               \
                                 return endian::big_to_native(*PTR_CAST(buffer, T))

/**
 * Converts a string (byte array) into a uint8_t.
 * @param str the string to convert
 * @return uint8_t value
 */
static uint8_t generate_uint8_of_string(const std::string &str) {
    INTEGRAL_OF_STRING(str, uint8_t);
}

/**
 * Converts a string (byte array) into a uint16_t.
 * @param str the string to convert
 * @return uint16_t value
 */
static uint16_t generate_uint16_of_string(const std::string &str) {
    INTEGRAL_OF_STRING(str, uint16_t);
}

/**
 * Converts a string (byte array) into a uint32_t.
 * @param str the string to convert
 * @return uint32_t value
 */
static uint32_t generate_uint32_of_string(const std::string &str) {
    INTEGRAL_OF_STRING(str, uint32_t);
}

/**
 * Read uint8_t from the connection.
 * @param server_connection the connection
 * @return read uint8_t
 */
static SerializableUINT8 read_u8(Server::Connection &server_connection) {
    auto message = server_connection.read(U8_SIZE);
    return SerializableUINT8{generate_uint8_of_string(message)};
}

/**
 * Read uint16_t from the connection.
 * @param server_connection the connection
 * @return read uint16_t
 */
static SerializableUINT16 read_u16(Server::Connection &server_connection) {
    auto message = server_connection.read(U16_SIZE);
    return SerializableUINT16{generate_uint16_of_string(message)};
}

/**
 * Read uint32_t from the connection.
 * @param server_connection the connection
 * @return read uint32_t
 */
static SerializableUINT32 read_u32(Server::Connection &server_connection) {
    auto message = server_connection.read(U32_SIZE);
    return SerializableUINT32{generate_uint32_of_string(message)};
}

/**
 * Read PlayerId from the connection.
 * @param server_connection the connection
 * @return read PlayerId
 */
static auto read_player_id = read_u8;

/**
 * Read String from the connection.
 * @param server_connection the connection
 * @return read String
 */
static SerializableString read_string(Server::Connection &server_connection) {
    auto len = read_u8(server_connection);
    auto str = server_connection.read(len._number);
    return SerializableString{str};
}

/**
 * Read Player from the connection.
 * @param server_connection the connection
 * @return read Player
 */
static SerializablePlayer read_player(Server::Connection &server_connection) {
    auto name = read_string(server_connection);
    auto address = read_string(server_connection);
    return { name, address };
}

/**
 * Read BombId from the connection.
 * @param server_connection the connection
 * @return read BombId
 */
static auto read_bomb_id = read_u32;

/**
 * Read Position from the connection.
 * @param server_connection the connection
 * @return read Position
 */
static SerializablePosition read_position(Server::Connection &server_connection) {
    auto x = read_u16(server_connection);
    auto y = read_u16(server_connection);
    return { x, y };
}

/**
 * Read Bomb from the connection.
 * @param server_connection the connection
 * @return read Bomb
 */
[[maybe_unused]] static SerializableBomb read_bomb(Server::Connection &server_connection) {
    auto position = read_position(server_connection);
    auto timer = read_u16(server_connection);
    return { position, timer };
}

/**
 * Read Score from the connection.
 * @param server_connection the connection
 * @return read Score
 */
static auto read_score = read_u32;

/**
 * Read SerializableList<T> from the connection.
 * @tparam T the type of a list
 * @param server_connection the connection
 * @param H the function reading the type T from connection
 * @return read SerializableList<T>
 */
template<Derived<Serializable> T>
static SerializableList<T> read_list(Server::Connection &server_connection, std::function<T(Server::Connection &)> H) {
    auto len = read_u32(server_connection);

    SerializableList<T> _list;
    for (size_t i = 0; i < len._number; ++i) {
        auto element = H(server_connection);
        _list.list.push_back(element);
    }

    return _list;
}

/**
 * Read List<T> from the connection.
 * @tparam T the type of a list
 * @param server_connection the connection
 * @param H the function reading the type T from connection
 * @return read List<T>
 */
template<typename T>
static List<T> read_simple_list(Server::Connection &server_connection, std::function<T(Server::Connection &)> H) {
    auto len = read_u32(server_connection)._number;

    List<T> _list;
    for (size_t i = 0; i < len; ++i) {
        auto element = H(server_connection);
        _list.push_back(element);
    }

    return _list;
}

/**
 * Read SerializableMap<T, U> from the connection.
 * @tparam T key type of a map
 * @tparam U value type of a map
 * @param server_connection the connection
 * @param H the function reading the type T from connection
 * @param P the function reading the type U from connection
 * @return Read SerializableMap<T, U>
 */
template<Derived<Serializable> T, Derived<Serializable> U>
static SerializableMap<T, U> read_map(Server::Connection &server_connection, std::function<T(Server::Connection &)> H,
                   std::function<U(Server::Connection &)> P) {

    auto len = read_u32(server_connection);

    SerializableMap<T, U> _map;
    for (size_t i = 0; i < len._number; ++i) {
        auto key = H(server_connection);
        auto value = P(server_connection);
        _map.map[key] = value;
    }

    return _map;
}

/**
 * Read BombPlaced from the connection.
 * @param server_connection the connection
 * @return read BombPlaced
 */
static EventPtr read_bomb_placed(Server::Connection &server_connection) {
    auto id = read_bomb_id(server_connection);
    auto position = read_position(server_connection);
    return std::make_shared<BombPlaced>(id, position);
}

/**
 * Read BombExploded from the connection.
 * @param server_connection the connection
 * @return read BombExploded
 */
static EventPtr read_bomb_exploded(Server::Connection &server_connection) {
    auto id = read_bomb_id(server_connection);
    auto robots_destroyed = read_list<SerializablePlayerId>(
            server_connection, read_player_id);
    auto blocks_destroyed = read_list<SerializablePosition>(
            server_connection, read_position);
    return std::make_shared<BombExploded>(id, robots_destroyed, blocks_destroyed);
}

/**
 * Read PlayerMoved from the connection.
 * @param server_connection the connection
 * @return read PlayerMoved
 */
static EventPtr read_player_moved(Server::Connection &server_connection) {
    auto id = read_player_id(server_connection);
    auto position = read_position(server_connection);
    return std::make_shared<PlayerMoved>(id, position);
}

/**
 * Read BlockPlaced from the connection.
 * @param server_connection the connection
 * @return read BlockPlaced
 */
static EventPtr read_block_placed(Server::Connection &server_connection) {
    auto position = read_position(server_connection);
    return std::make_shared<BlockPlaced>(position);
}

using EventHandler = std::function<EventPtr(Server::Connection &)>;

/**
 * Read Event from the connection.
 * @param server_connection the connection
 * @return read Event
 */
static EventPtr read_event(Server::Connection &server_connection) {
    static EventHandler event_handlers[] = {
            read_bomb_placed,   // 0
            read_bomb_exploded, // 1
            read_player_moved,  // 2
            read_block_placed   // 3
    };
    static size_t handlers = sizeof(event_handlers) / sizeof(event_handlers[0]);

    auto type = read_u8(server_connection);
    if (type._number >= handlers)
        throw DeserializationException{};

    return event_handlers[type._number](server_connection);
}

/**
 * Read Hello Message from the connection.
 * @param server_connection the connection
 * @return read Hello Message
 */
static ServerMessagePtr read_hello(Server::Connection &server_connection) {
    auto server_name = read_string(server_connection);
    auto players_count = read_u8(server_connection);
    auto size_x = read_u16(server_connection);
    auto size_y = read_u16(server_connection);
    auto game_length = read_u16(server_connection);
    auto explosion_radius = read_u16(server_connection);
    auto bomb_timer = read_u16(server_connection);

    return std::make_shared<Hello>(
            server_name,
            players_count,
            size_x, size_y,
            game_length,
            explosion_radius,
            bomb_timer);
}

/**
 * Read AcceptedPlayer Message from the connection.
 * @param server_connection the connection
 * @return read AcceptedPlayer Message
 */
static ServerMessagePtr read_accepted_player(Server::Connection &server_connection) {
    auto id = read_player_id(server_connection);
    auto player = read_player(server_connection);
    return std::make_shared<AcceptedPlayer>(id, player);
}

/**
 * Read GameStarted Message from the connection.
 * @param server_connection the connection
 * @return read GameStarted Message
 */
static ServerMessagePtr read_game_started(Server::Connection &server_connection) {
    auto players = read_map<SerializablePlayerId, SerializablePlayer>(
            server_connection,
            read_player_id,
            read_player);
    return std::make_shared<GameStarted>(players);
}

/**
 * Read Turn Message from the connection.
 * @param server_connection the connection
 * @return read Turn Message
 */
static ServerMessagePtr read_turn(Server::Connection &server_connection) {
    auto turn = read_u16(server_connection);
    auto events = read_simple_list<EventPtr>(server_connection, read_event);
    return std::make_shared<Turn>(turn, events);
}

/**
 * Read GameEnded Message from the connection.
 * @param server_connection the connection
 * @return read GameEnded Message
 */
static ServerMessagePtr read_game_ended(Server::Connection &server_connection) {
    auto scores = read_map<SerializablePlayerId, SerializableScore>(
            server_connection,
            read_player_id,
            read_score);
    return std::make_shared<GameEnded>(scores);
}

using ServerMessageHandler = std::function<ServerMessagePtr(Server::Connection &)>;

ServerMessagePtr read_server_message(Server::Connection &server_connection) {
    static ServerMessageHandler server_message_handlers[] = {
            read_hello,             // 0
            read_accepted_player,   // 1
            read_game_started,      // 2
            read_turn,              // 3
            read_game_ended         // 4
    };

    auto message_type = read_u8(server_connection)._number;

    size_t handlers = sizeof(server_message_handlers) / sizeof(server_message_handlers[0]);
    if (message_type >= handlers)
        throw DeserializationException{};

    return server_message_handlers[message_type](server_connection);
}
