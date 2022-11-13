#ifndef SIK_BOMBERMAN_SERVERMESSAGE_H
#define SIK_BOMBERMAN_SERVERMESSAGE_H

#include <utility>
#include <iostream>
#include <iomanip>

#include "Networkable.h"
#include "String.h"
#include "Numbers.h"
#include "Player.h"
#include "Containers.h"
#include "Bomb.h"

namespace Models {

    class ServerMessage : public Networkable {
    public:
        void deserialize(Client::Connection &connection) override {
            (void) connection;
            assert(false); // Should never be used.
        }
    };

    class Hello : public ServerMessage {
        constexpr static const UINT8 ID = 0;
    public:
        Hello() = default;

        Hello(String serverName,
              UINT8 playersCount,
              UINT16 sizeX, UINT16 sizeY,
              UINT16 gameLength,
              UINT16 explosionRadius,
              UINT16 bombTimer) :
                  server_name(std::move(serverName)),
                  players_count(std::move(playersCount)),
                  size_x(std::move(sizeX)), size_y(std::move(sizeY)),
                  game_length(std::move(gameLength)),
                  explosion_radius(std::move(explosionRadius)),
                  bomb_timer(std::move(bombTimer)) {}

        [[nodiscard]] std::string serialize() const override {
            return ID.serialize() +
                server_name.serialize() +
                players_count.serialize() +
                size_x.serialize() +
                size_y.serialize() +
                game_length.serialize() +
                explosion_radius.serialize() +
                bomb_timer.serialize();
        }

        String server_name;
        UINT8 players_count;
        UINT16 size_x;
        UINT16 size_y;
        UINT16 game_length;
        UINT16 explosion_radius;
        UINT16 bomb_timer;
    };

    class AcceptedPlayer : public ServerMessage {
        constexpr static const UINT8 ID = 1;
    public:
        AcceptedPlayer() = default;

        AcceptedPlayer(PlayerId id, Player player) : id(std::move(id)), player(std::move(player)) {}

        [[nodiscard]] std::string serialize() const override {
            return ID.serialize() + id.serialize() + player.serialize();
        }

        PlayerId id;
        Player player;
    };

    class GameStarted : public ServerMessage {
        constexpr static const UINT8 ID = 2;
    public:
        GameStarted() = default;

        explicit GameStarted(const Map<PlayerId, Player> &players) : players(players) {}

        [[nodiscard]] std::string serialize() const override {
            return ID.serialize() + players.serialize();
        }

        Map<PlayerId, Player> players;
    };

    class Event : public Networkable {
    public:
        void deserialize(Client::Connection &connection) override {
            (void) connection;
            assert(false); // Should never be used.
        }
    };

    class BombPlaced : public Event {
        constexpr static const UINT8 ID = 0;
    public:
        BombPlaced() = default;

        BombPlaced(BombId id, Position position) : id(std::move(id)), position(std::move(position)) {}

        [[nodiscard]] std::string serialize() const override {
            LOG(debug) << "Serializing BombPlaced id=" << id.value << " position=" << position.x.value << "," << position.y.value;
            return ID.serialize() + id.serialize() + position.serialize();
        }

        BombId id;
        Position position;
    };

    class BombExploded : public Event {
        constexpr static const UINT8 ID = 1;
    public:
        [[nodiscard]] std::string serialize() const override {
            return ID.serialize() +
                id.serialize() +
                robots_destroyed.serialize() +
                blocks_destroyed.serialize();
        }

        BombId id;
        List<PlayerId> robots_destroyed;
        List<Position> blocks_destroyed;
    };

    class PlayerMoved : public Event {
        constexpr static const UINT8 ID = 2;
    public:
        PlayerMoved() = default;

        PlayerMoved(PlayerId id, Position position) : id(std::move(id)), position(std::move(position)) {}

        [[nodiscard]] std::string serialize() const override {
            return ID.serialize() + id.serialize() + position.serialize();
        }

        PlayerId id;
        Position position;
    };

    class BlockPlaced : public Event {
        constexpr static const UINT8 ID = 3;
    public:
        BlockPlaced() = default;

        BlockPlaced(Position position) : position(std::move(position)) {}

        [[nodiscard]] std::string serialize() const override {
            return ID.serialize() + position.serialize();
        }

        Position position;
    };

    class Turn : public ServerMessage {
        constexpr static const UINT8 ID = 3;
    public:
        [[nodiscard]] std::string serialize() const override {
            return ID.serialize() + turn.serialize() + events.serialize();
        }

        UINT16 turn;
        PtrList<Event> events;
    };

    class GameEnded : public ServerMessage {
        constexpr static const UINT8 ID = 4;
    public:
        GameEnded() {}

        GameEnded(const Map <PlayerId, Score> &scores) : scores(scores) {}

        std::string serialize() const override {
            return ID.serialize() + scores.serialize();
        }

        Map<PlayerId, Score> scores;
    };
}


#endif //SIK_BOMBERMAN_SERVERMESSAGE_H
