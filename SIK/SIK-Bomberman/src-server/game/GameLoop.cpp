#include "GameLoop.h"
#include "models/models.h"

#include <utility>

void Game::GameLoop::lobby() {
    game_state.is_lobby = true;
    std::vector<std::shared_ptr<ServerMessage>> messages;
    instances->for_each([&](std::shared_ptr<PlayerInstance> &instance) {
        instance->input_actions.perform([&](std::queue<std::shared_ptr<ClientMessage>> &queue) {
            while (!queue.empty()) {
                auto message = queue.front();
                queue.pop();
                if (message->get_id() != JOIN_ID || instance->is_playing) { continue; }
                if (game_state.players.map.size() != *arguments->players_count) {
                    message->update_state(game_state, PlayerId{}/* unused anyway */);
                    auto new_player = *game_state.players.map.rbegin();
                    auto player = new_player.second;
                    auto id = new_player.first;
                    instance->is_playing = true;
                    instance->player_id = id;

                    auto player_accepted = std::make_shared<Models::AcceptedPlayer>(id, player);
                    messages.push_back(player_accepted);
                }
            }
        });
    });

    for (auto &message : messages) {
        instances->for_each([&](std::shared_ptr<PlayerInstance> &instance) {
            instance->output_actions.push(message);
        });
    }
}

/**
 * Calculate x mod. m for x being uint32 and m being uint16.
 * @param x uint32
 * @param m uint16
 * @return x mod. m
 */
static uint16_t mod_32_16(uint32_t x, uint16_t m) {
    return static_cast<uint16_t>(x % ((uint32_t) m));
}

void Game::GameLoop::init_game() {
    auto turn = std::make_shared<Models::Turn>();
    turn->turn = 0;
    auto &events = turn->events.list;

    for (auto &[id, player] : game_state.players.map) {
        player.position.x = mod_32_16(rng->value(), *arguments->size_x);
        player.position.y = mod_32_16(rng->value(), *arguments->size_y);
        events.push_back(std::make_shared<Models::PlayerMoved>(id, player.position));
    }

    for (uint16_t i = 0; i < *arguments->initial_blocks; ++i) {
        Position pos;
        pos.x = mod_32_16(rng->value(), *arguments->size_x);
        pos.y = mod_32_16(rng->value(), *arguments->size_y);
        if (game_state.blocks.insert(pos).second) {
            events.push_back(std::make_shared<Models::BlockPlaced>(pos));
        }
    }

    instances->for_each([&](std::shared_ptr<PlayerInstance> &instance) {
        instance->output_actions.push(turn);
    });
    game_state.current_messages.push(turn);
}

void Game::GameLoop::calculate_destruction(const Models::Position& position,
                                           const std::shared_ptr<BombExploded>& bomb_exploded) {

    static int dx[] = { 0, 1, 0, -1 };
    static int dy[] = { 1, 0, -1, 0 };
    static size_t ds = 4;

    std::set<Models::Position> explosion;

    for (size_t k = 0; k < ds; ++k) {
        int current_x = position.x.value;
        int current_y = position.y.value;

        for (size_t e = 0; e <= *arguments->explosion_radius; ++e) {
            if (current_x < 0 || current_y < 0 || current_x >= *arguments->size_x || current_y >= *arguments->size_y) {
                break;
            }

            auto current = Models::Position{current_x, current_y};
            explosion.insert(current);
            auto block_it = game_state.blocks.find(current);
            if (game_state.blocks.end() != block_it) {
                game_state.blocks_to_destroy.insert(current);
                bomb_exploded->blocks_destroyed.list.push_back(current);
                break;
            }

            current_x += dx[k];
            current_y += dy[k];
        }
    }

    // Remove duplicates from blocks_destroyed
    std::set<Models::Position> _blocks{ bomb_exploded->blocks_destroyed.list.begin(), bomb_exploded->blocks_destroyed.list.end() };
    bomb_exploded->blocks_destroyed.list.clear();
    bomb_exploded->blocks_destroyed.list.insert(bomb_exploded->blocks_destroyed.list.begin(), _blocks.begin(), _blocks.end());

    for (auto &[player_id, player] : game_state.players.map) {
        if (explosion.contains(player.position)) {
            game_state.players_to_destroy.insert(player_id);
            bomb_exploded->robots_destroyed.list.push_back(player_id);
        }
    }
}

std::shared_ptr<Models::ServerMessage> Game::GameLoop::update_game(uint16_t turn) {
    std::map<PlayerId, std::shared_ptr<Models::ClientMessage>> last_reasonable_move;

    instances->for_each([&](std::shared_ptr<PlayerInstance> &instance) {
        while (!instance->input_actions.empty()) {
            auto message = instance->input_actions.take();
            if (message->get_id() != JOIN_ID) {
                last_reasonable_move[instance->player_id] = message;
            }
        }
    });

    std::shared_ptr<Models::Turn> turn_message = std::make_shared<Models::Turn>();
    turn_message->turn = turn;
    auto &events = turn_message->events;

    std::vector<BombId> bombs_to_remove;
    for (auto &[bomb_id, bomb] : game_state.bombs) {
        auto diff_zero = --bomb.timer <=> 0;
        if (diff_zero == 0) {
            auto bomb_explosion = std::make_shared<Models::BombExploded>();
            bomb_explosion->id = bomb_id;
            bombs_to_remove.push_back(bomb_id);
            calculate_destruction(bomb.position, bomb_explosion);
            events.list.push_back(bomb_explosion);
        }
    }
    for (auto &bomb_id : bombs_to_remove) {
        game_state.bombs.erase(bomb_id);
    }

    for (auto &block : game_state.blocks_to_destroy) {
        game_state.blocks.erase(block);
    }

    for (auto &[player_id, player] : game_state.players.map) {
        game_state.scores.map[player_id];
        if (game_state.players_to_destroy.contains(player_id)) {
            // zniszczony
            player.position.x = mod_32_16(rng->value(), *arguments->size_x);
            player.position.y = mod_32_16(rng->value(), *arguments->size_y);
            events.list.push_back(std::make_shared<Models::PlayerMoved>(player_id, player.position));
            ++game_state.scores.map[player_id];
        }
        else {
            // nie zniszczony
            auto player_move = last_reasonable_move.find(player_id);
            if (player_move != last_reasonable_move.end()) {
                // move was played
                if (player_move->second->get_id() != JOIN_ID) {
                    auto event = player_move->second->update_state(game_state, player_id);
                    if (event) {
                        events.list.push_back(event);
                    }
                }
            }
        }
    }

    game_state.blocks_to_destroy.clear();
    game_state.players_to_destroy.clear();

    return turn_message;
}

void Game::GameLoop::game() {
    game_state.is_lobby = false;

    init_game();
    // Send turn=0

    auto round_time = std::chrono::milliseconds(*arguments->turn_duration);
    for (uint16_t turn = 1; turn <= *arguments->game_length; ++turn) {
        auto before = std::chrono::system_clock::now();
        auto turn_message = update_game(turn);
        instances->for_each([&](std::shared_ptr<PlayerInstance> &instance) {
            instance->input_actions.clear();
            instance->output_actions.push(turn_message);
        });
        game_state.current_messages.push(turn_message);
        auto after = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(after - before);
        if (turn != *arguments->game_length)
            std::this_thread::sleep_for(round_time - elapsed);
    }
}

void Game::GameLoop::clear_state() {
    game_state.players.map.clear();
    game_state.bombs.clear();
    game_state.blocks.clear();
    game_state.scores.map.clear();

    game_state.blocks_to_destroy.clear();
    game_state.players_to_destroy.clear();

    game_state.bomb_counter = 0;
    game_state.current_messages.clear();

    instances->for_each([&](std::shared_ptr<PlayerInstance> &instance) {
        instance->input_actions.clear();
        instance->is_playing = false;
    });
}

void Game::GameLoop::operator()() {
    std::mutex mutex;
    std::condition_variable wait_for_input;
    bool wake{false};

    instances->set_listener([&] {
        wake = true;
        wait_for_input.notify_all();
    });

    try_forever(
        while (true) {
            lobby();
            if (game_state.players.map.size() != *arguments->players_count) {
                std::unique_lock<std::mutex> lock{mutex};
                wait_for_input.wait_for(lock, std::chrono::seconds(1), [&]() {
                    return !wake;
                });
                wake = false;
            }
            else {
                break;
            }
        }

        auto game_started = std::make_shared<Models::GameStarted>(game_state.players);
        instances->for_each([&](std::shared_ptr<PlayerInstance> &instance) {
            instance->output_actions.push(game_started);
        });
        game_state.current_messages.push(game_started);

        LOG(debug) << "Starting game";
        game();
        LOG(debug) << "Game finished";

        auto game_ended = std::make_shared<Models::GameEnded>(game_state.scores);
        instances->for_each([&](std::shared_ptr<PlayerInstance> &instance) {
            instance->input_actions.clear();
            instance->output_actions.push(game_ended);
        });

        instances->lock();
        std::vector<std::shared_ptr<PlayerInstance>> to_remove;
        instances->for_each_unsafe([&](std::shared_ptr<PlayerInstance> &instance) {
            LOG(debug) << "is disconnected: " << instance->disconnected;
            if (instance->disconnected) {
                to_remove.push_back(instance);
            }
        });
        instances->perform([&](std::set<std::shared_ptr<PlayerInstance>> &_instances) {
            for (auto &instance : to_remove) {
                connections.remove(receivers[instance]->connection);
                senders.erase(instance);
                receivers.erase(instance);
                _instances.erase(instance);
            }
        });

        instances->unlock();

        clear_state();
    )
    catch (const std::exception &exception) {
        LOG(debug) << "Exception in the game. " << exception.what();
    }
}

Game::GameLoop::GameLoop(std::shared_ptr<Server::Arguments> _arguments,
                         std::shared_ptr<safe::set<std::shared_ptr<PlayerInstance>>> instances,
                         std::shared_ptr<Server::RNG> rng,
                         safe::queue<std::shared_ptr<Models::ServerMessage>> &current_messages,
                         std::map<std::shared_ptr<PlayerInstance>, std::shared_ptr<Client::Sender>> &senders,
                         std::map<std::shared_ptr<PlayerInstance>, std::shared_ptr<Client::Receiver>> &receivers,
                         safe::set<std::shared_ptr<Client::Connection>> &connections)
        : arguments(std::move(_arguments)), instances(std::move(instances)),
        rng(std::move(rng)), game_state(arguments, current_messages),
        senders(senders), receivers(receivers), connections(connections) {}

