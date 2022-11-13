#ifndef SIK_BOMBERMAN_CLIENTMESSAGE_H
#define SIK_BOMBERMAN_CLIENTMESSAGE_H


#include "Networkable.h"
#include "String.h"
#include "Direction.h"
#include "../GameState.h"
#include "../../client/Connection.h"
#include "../../utils/utils.h"
#include "ServerMessage.h"

namespace Models {

    class ClientMessage : public Networkable {
    public:
        explicit ClientMessage(uint8_t id) : id(id) {}
        virtual std::shared_ptr<Models::Event> update_state(Game::GameState &state,
                                                                    Models::PlayerId player_id) = 0;

        [[nodiscard]] std::string serialize() const override {
            assert(false);
            return ""; // never used
        }

        [[nodiscard]] uint8_t get_id() const {
            return id;
        }

    private:
        uint8_t id;
    };

    class Join : public ClientMessage {
    public:
        explicit Join(uint8_t id) : ClientMessage(id) {}

        void deserialize(Client::Connection &connection) override {
            name.deserialize(connection);
            address.value = connection.get_endpoint();
        }

        std::shared_ptr<Models::Event> update_state(Game::GameState &state, Models::PlayerId player_id) override {
            LOG(debug) << "Updating Join.";
            (void) player_id;

            auto new_id = static_cast<uint8_t>(state.players.map.size());
            auto player = Player{ name, address };
            state.players.map[new_id] = player;
            state.scores.map[new_id] = Score{0};
            return nullptr;
        }

        String name;
        String address;
    };

    class PlaceBomb : public ClientMessage {
    public:
        explicit PlaceBomb(uint8_t id) : ClientMessage(id) {}

        void deserialize(Client::Connection &) override {}

        std::shared_ptr<Models::Event> update_state(Game::GameState &state, Models::PlayerId player_id) override {
            LOG(debug) << "Updating PlaceBomb.";
            (void) state;

            auto position = state.players.map[player_id].position;
            auto new_id = BombId{++state.bomb_counter};
            state.bombs[new_id] = Models::Bomb{ position, *state.arguments->bomb_timer };
            return std::make_shared<Models::BombPlaced>(new_id, position);
        }
    };

    class PlaceBlock : public ClientMessage {
    public:
        explicit PlaceBlock(uint8_t id) : ClientMessage(id) {}

        void deserialize(Client::Connection &) override {}

        std::shared_ptr<Models::Event> update_state(Game::GameState &state, Models::PlayerId player_id) override {
            LOG(debug) << "Updating BlockPlaced.";

            auto position = state.players.map[player_id].position;
            if (!state.blocks.contains(position)) {
                state.blocks.insert(position);
                return std::make_shared<Models::BlockPlaced>(position);
            }
            else {
                return nullptr;
            }
        }
    };

    class Move : public ClientMessage {
    public:
        explicit Move(uint8_t id) : ClientMessage(id) {}

        void deserialize(Client::Connection &connection) override {
            direction.deserialize(connection);
        }

        std::shared_ptr<Models::Event> update_state(Game::GameState &state, Models::PlayerId player_id) override {
            LOG(debug) << "Updating Move.";
            auto position = state.players.map[player_id].position;
            bool overflow = false;
            switch (direction.dir) {
                case Direction::Up:
                    if (position.y.value + 1 == *state.arguments->size_y) overflow = true;
                    position.y += 1;
                    break;
                case Direction::Right:
                    if (position.x.value + 1 == *state.arguments->size_x) overflow = true;
                    position.x += 1;
                    break;
                case Direction::Down:
                    if (position.y.value == 0) overflow = true;
                    position.y -= 1;
                    break;
                case Direction::Left:
                    if (position.x.value == 0) overflow = true;
                    position.x -= 1;
                    break;
            }
            if (!overflow && !state.blocks.contains(position)) {
                state.players.map[player_id].position = position;
                return std::make_shared<Models::PlayerMoved>(player_id, position);
            }
            else {
                return nullptr;
            }
        }

        Direction direction;
    };

    template<std::derived_from<ClientMessage> T>
    static auto handle(uint8_t id) {
        return [id](Client::Connection &connection) {
            std::shared_ptr<ClientMessage> message = std::make_shared<T>(id);
            message->deserialize(connection);
            return message;
        };
    }

    [[maybe_unused]] static const uint8_t JOIN_ID = 0;
    [[maybe_unused]] static const uint8_t PLACE_BOMB_ID = 1;
    [[maybe_unused]] static const uint8_t PLACE_BLOCK_ID = 2;
    [[maybe_unused]] static const uint8_t MOVE_ID = 3;

    [[maybe_unused]] static std::shared_ptr<ClientMessage> ClientMessageFactory(Client::Connection &connection) {
        UINT8 message_id_;
        message_id_.deserialize(connection);
        size_t message_id = message_id_.value;

        static std::function<std::shared_ptr<ClientMessage>(Client::Connection &)> handlers[] = {
                handle<Join>(JOIN_ID),
                handle<PlaceBomb>(PLACE_BOMB_ID),
                handle<PlaceBlock>(PLACE_BLOCK_ID),
                handle<Move>(MOVE_ID)
        };
        static size_t handlers_count = sizeof(handlers) / sizeof(handlers[0]);

        if (message_id >= handlers_count)
            throw DeserializationException{};

        return handlers[message_id](connection);
    }
}


#endif //SIK_BOMBERMAN_CLIENTMESSAGE_H
