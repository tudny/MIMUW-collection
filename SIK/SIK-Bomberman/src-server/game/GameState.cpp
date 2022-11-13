#include "GameState.h"

#include <utility>

Game::GameState::GameState(std::shared_ptr<Server::Arguments> arguments,
                           safe::queue<std::shared_ptr<Models::ServerMessage>> &current_messages)
                           : arguments(std::move(arguments)), current_messages(current_messages) {}
