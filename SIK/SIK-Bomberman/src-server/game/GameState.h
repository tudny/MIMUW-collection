#ifndef SIK_BOMBERMAN_GAMESTATE_H
#define SIK_BOMBERMAN_GAMESTATE_H


#include <set>
#include "models/Player.h"
#include "models/Position.h"
#include "models/Bomb.h"
#include "models/Containers.h"
#include "models/ServerMessage.h"
#include "../utils/arguments.h"
#include "../threads/queue.h"

namespace Game {

    using namespace Models;

    class GameState {
    public:
        explicit GameState(std::shared_ptr<Server::Arguments> arguments, safe::queue<std::shared_ptr<Models::ServerMessage>> &current_messages);

        bool is_lobby{true};
        Map<PlayerId, Player> players;
        std::map<BombId, Bomb> bombs;
        uint32_t bomb_counter = 0;
        std::set<Position> blocks;
        Map<PlayerId, Score> scores;

        std::set<Position> blocks_to_destroy;
        std::set<PlayerId> players_to_destroy;

        std::shared_ptr<Server::Arguments> arguments;

        safe::queue<std::shared_ptr<Models::ServerMessage>> &current_messages;
    };
}


#endif //SIK_BOMBERMAN_GAMESTATE_H
