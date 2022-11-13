#ifndef SIK_BOMBERMAN_GAMELOOP_H
#define SIK_BOMBERMAN_GAMELOOP_H

#include <memory>
#include <set>
#include "../utils/arguments.h"
#include "models/models.h"
#include "../threads/set.h"
#include "../client/PlayerInstance.h"
#include "../utils/RNG.h"
#include "../client/Sender.h"
#include "../client/Receiver.h"

namespace Game {

    using namespace Models;

    class GameLoop {
    public:
        GameLoop(std::shared_ptr<Server::Arguments> arguments,
                 std::shared_ptr<safe::set<std::shared_ptr<PlayerInstance>>> instances,
                 std::shared_ptr<Server::RNG> rng,
                 safe::queue<std::shared_ptr<Models::ServerMessage>> &current_messages,
                 std::map<std::shared_ptr<PlayerInstance>, std::shared_ptr<Client::Sender>> &senders,
                 std::map<std::shared_ptr<PlayerInstance>, std::shared_ptr<Client::Receiver>> &receivers,
                 safe::set<std::shared_ptr<Client::Connection>> &connections);

        void operator()();

    private:
        void lobby();
        void game();

        void init_game();
        std::shared_ptr<Models::ServerMessage> update_game(uint16_t turn);
        void clear_state();

        std::shared_ptr<Server::Arguments> arguments;
        std::shared_ptr<safe::set<std::shared_ptr<PlayerInstance>>> instances;
        std::shared_ptr<Server::RNG> rng;
        GameState game_state;

        void calculate_destruction(const Models::Position& position, const std::shared_ptr<BombExploded>& bomb_exploded);

        std::map<std::shared_ptr<PlayerInstance>, std::shared_ptr<Client::Sender>> &senders;
        std::map<std::shared_ptr<PlayerInstance>, std::shared_ptr<Client::Receiver>> &receivers;
        safe::set<std::shared_ptr<Client::Connection>> &connections;
    };
}


#endif //SIK_BOMBERMAN_GAMELOOP_H
