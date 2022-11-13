#ifndef SIK_BOMBERMAN_PLAYERINSTANCE_H
#define SIK_BOMBERMAN_PLAYERINSTANCE_H


#include <memory>
#include "../game/models/models.h"
#include "../threads/queue.h"

class PlayerInstance {
public:
    bool is_playing{false};
    bool disconnected{false};
    Models::PlayerId player_id;
    safe::queue<std::shared_ptr<Models::ClientMessage>> input_actions;
    safe::queue<std::shared_ptr<Models::ServerMessage>> output_actions;
};


#endif //SIK_BOMBERMAN_PLAYERINSTANCE_H
