#ifndef SIK_BOMBERMAN_GAMESTATE_H
#define SIK_BOMBERMAN_GAMESTATE_H

/* STD */
#include <string>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <cinttypes>

/* LOCAL */
#include "../serializer/serializer.h"
#include "../server/Connection.h"

/* SIMPLE TYPES */

template<typename T>
using List = std::vector<T>;

/**
 * Represents the state od the game.
 */
struct GameState {
    /**
     * Default constructor. Sets all values to default values.
     */
    GameState();

    // Server info
    SerializableString server_name;
    SerializableUINT8 players_count{};
    SerializableUINT16 size_x{};
    SerializableUINT16 size_y{};
    SerializableUINT16 game_length{};
    SerializableUINT16 explosion_radius{};
    SerializableUINT16 bomb_timer{};

    // Server status
    bool is_lobby{true};
    bool do_send{true};

    // Data
    SerializableUINT16 turn{};
    SerializableMap<SerializablePlayerId, SerializablePlayer> players;
    SerializableMap<SerializablePlayerId, SerializablePosition> player_position;
    std::set<SerializablePosition> blocks;
    SerializableMap<SerializableBombId, SerializableBomb> bombs; // when sending use as List<Bomb>
    std::set<SerializablePosition> explosions;
    SerializableMap<SerializablePlayerId, SerializableScore> scores;

    // Temporary date
    std::set<SerializablePlayerId> players_destroyed_in_turn;
    std::set<SerializablePosition> blocks_to_remove;
};

/**
 * Converts a set of type <T> into a SerializableList<T>.
 * @tparam T the type of a set. It must derive from Serializable.
 * @param _set set to be converted
 * @return SerializableList<T> with the elements from original set
 */
template<Derived<Serializable> T>
SerializableList<T> from_set(std::set<T> _set) {
    return SerializableList{std::vector<T>{_set.begin(), _set.end()}};
}

/**
 * Converts a map of type <U, T> into a SerializableList<T> from values
 * set of the map.
 * @tparam U map key type
 * @tparam T map value type
 * @param _map map of which values should be converted into a list
 * @return a list of values of the map as SerializableList<T>
 */
template<Derived<Serializable> U, Derived<Serializable> T>
SerializableList<T> from_map_values(SerializableMap<U, T> _map) {
    std::vector<T> _list;
    for (auto &[_, value]: _map.map) {
        _list.push_back(value);
    }
    return SerializableList<T>{_list};
}


/* EVENTS */

/**
 * An abstract class representing an Event.
 */
class Event {
    /**
     * An abstract method called on an instance of an Event
     * to change the game state accordingly.
     * @param game_state the state of game to change
     */
    virtual void e_handler(GameState &game_state) = 0;

public:
    /**
     * A method called to change the state of the game.
     * It does some work predefined for each event,
     * so inside it callas, so far, abstract e_handler.
     * @param game_state the state of game to change.
     */
    void handle(GameState &game_state);
};

/**
 * Class representing the event of placing a bomb.
 */
class BombPlaced : public Event {
    /**
     * An implementation of bomb placing event.
     * @param game_state the state of game to change
     */
    void e_handler(GameState &game_state) override;

public:
    /**
     * A constructor of BombPlaced event.
     * @param id id of a bomb
     * @param position a position of placed bomb
     */
    BombPlaced(const SerializableBombId &id, const SerializablePosition &position);

    SerializableBombId id; ///< id of a bomb
    SerializablePosition position; ///< position of a placed bomb
};

/**
 * Class representing the event of the bomb explosion.
 */
class BombExploded : public Event {
    /**
     * An implementation of the bomb explosion.
     * @param game_state the state of game to change
     */
    void e_handler(GameState &game_state) override;

public:
    /**
     * A constructor of BombExploded event.
     * @param id is of a bomb
     * @param robotsDestroyed a list of destroyed robots
     * @param blocksDestroyed a list of destroyed blocks
     */
    BombExploded(const SerializableBombId &id,
                 const SerializableList<SerializablePlayerId> &robotsDestroyed,
                 const SerializableList<SerializablePosition> &blocksDestroyed);

    SerializableBombId id; ///< id of a bomb
    SerializableList<SerializablePlayerId> robots_destroyed; ///< list of destroyed robots
    SerializableList<SerializablePosition> blocks_destroyed; ///< list of destroyed blocks
};

/**
 * Class representing the event of player's move.
 */
class PlayerMoved : public Event {
    /**
     * An implementation of a player's move.
     * @param game_state the state of game to change
     */
    void e_handler(GameState &game_state) override;

public:
    /**
     * A constructor of PlayerMoved event.
     * @param id id of a player
     * @param position new position of a plyer
     */
    PlayerMoved(const SerializablePlayerId &id, const SerializablePosition &position);

    SerializablePlayerId id; ///< id of a player
    SerializablePosition position; ///< new position of a player
};

/**
 * Class representing the event of BlockPlaced.
 */
class BlockPlaced : public Event {
    /**
     * An implementation of a placed block.
     * @param game_state the state of game to change
     */
    void e_handler(GameState &game_state) override;

public:
    /**
     * A constructor of BlockPlaced event.
     * @param position a position of a new block
     */
    explicit BlockPlaced(const SerializablePosition &position);

    SerializablePosition position; ///< a position of a new block
};

using EventPtr = std::shared_ptr<Event>;


/* SERVER MESSAGE */

/**
 * A class representing a message from the server.
 */
class ServerMessage {
public:
    /**
     * An abstract handler of a message. It changes the game state.
     * @param game_state the state of a game to be changed
     */
    virtual void handle(GameState &game_state) = 0;
};

/**
 * An implementation of a Hello Message.
 */
class Hello : public ServerMessage {
public:
    Hello(const SerializableString &serverName,
          const SerializableUINT8 &playersCount,
          const SerializableUINT16 &sizeX,
          const SerializableUINT16 &sizeY,
          const SerializableUINT16 &gameLength,
          const SerializableUINT16 &explosionRadius,
          const SerializableUINT16 &bombTimer);

    void handle(GameState &game_state) override;

    SerializableString server_name;
    SerializableUINT8 players_count;
    SerializableUINT16 size_x;
    SerializableUINT16 size_y;
    SerializableUINT16 game_length;
    SerializableUINT16 explosion_radius;
    SerializableUINT16 bomb_timer;
};

/**
 * An implementation of a AcceptedPlayer Message.
 */
class AcceptedPlayer : public ServerMessage {
public:
    AcceptedPlayer(const SerializablePlayerId &id, const SerializablePlayer &player);

    void handle(GameState &game_state) override;

    SerializablePlayerId id{};
    SerializablePlayer player{};
};

/**
 * An implementation of a GameStarted Message.
 */
class GameStarted : public ServerMessage {
public:
    explicit GameStarted(const SerializableMap<SerializablePlayerId, SerializablePlayer> &players);

    void handle(GameState &game_state) override;

    SerializableMap<SerializablePlayerId, SerializablePlayer> players;
};

/**
 * An implementation of a Turn Message.
 */
class Turn : public ServerMessage {
public:
    Turn(const SerializableUINT16 &turn, const List<EventPtr> &events);

    void handle(GameState &game_state) override;

    SerializableUINT16 turn{};
    List<EventPtr> events{};
};

/**
 * An implementation of a Turn Message.
 */
class GameEnded : public ServerMessage {
public:
    explicit GameEnded(const SerializableMap<SerializablePlayerId, SerializableScore> &scores);

    void handle(GameState &game_state) override;

    SerializableMap<SerializablePlayerId, SerializableScore> scores{};
};

using ServerMessagePtr = std::shared_ptr<ServerMessage>;

namespace {
    enum {
        U8_SIZE = sizeof(uint8_t),
        U16_SIZE = sizeof(uint16_t),
        U32_SIZE = sizeof(uint32_t),
    };
}

/**
 * Read the whole message from server and parse it accordingly.
 * This function needs the connection in order to read appropriate
 * amount of bytes based on message type, structure length etc.
 * @param server_connection the connection to the server
 * @return a pointer to message a structure with appropriate handler
 */
ServerMessagePtr read_server_message(Server::Connection &server_connection);


#endif //SIK_BOMBERMAN_GAMESTATE_H
