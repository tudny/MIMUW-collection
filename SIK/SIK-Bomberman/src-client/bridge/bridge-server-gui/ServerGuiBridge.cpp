/* BOOST */
#include <boost/format.hpp>

/* LOCAL */
#include "ServerGuiBridge.h"

/**
 * Pulls out needed information from the game state to
 * create Lobby message.
 * and converts them into lobby draw message.
 * @param game_state the game state
 * @return pointer to lobby structure ready to be sent to GUI
 */
static std::shared_ptr<DrawMessage> generate_lobby_message(GameState &game_state) {
    return std::make_shared<Lobby>(
        game_state.server_name,
        game_state.players_count,
        game_state.size_x,
        game_state.size_y,
        game_state.game_length,
        game_state.explosion_radius,
        game_state.bomb_timer,
        game_state.players
    );
}

/**
 * Pulls out needed information from the game state to
 * create Game message.
 * @param game_state the game state
 * @return pointer to game structure ready to be sent to GUI
 */
static std::shared_ptr<DrawMessage> generate_game_message(GameState &game_state) {
    return std::make_shared<Game>(
        game_state.server_name,
        game_state.size_x,
        game_state.size_y,
        game_state.game_length,
        game_state.turn,
        game_state.players,
        game_state.player_position,
        from_set(game_state.blocks),
        from_map_values(game_state.bombs),
        from_set(game_state.explosions),
        game_state.scores
    );
}

/* BRIDGE */

ServerGuiBridge::ServerGuiBridge(
        GUI::Connection &gui_connection,
        Server::Connection &server_connection,
        const callback_t &callback,
        GameState &game_state
) :
        Bridge(gui_connection,
               server_connection,
               callback,
               game_state) {}

void ServerGuiBridge::operator()() {

    try_forever (
        LOG(debug) << "[SERVER->] Waiting for message from server.";
        auto message_op = read_server_message(server_connection_);
        LOG(debug) << "[SERVER->] Received message from server.";

        message_op->handle(game_state);

        LOG(debug) << boost::format("[->GUI] Sending to GUI. {is_lobby=%d}") % game_state.is_lobby;
        if (game_state.do_send) {
            auto message_ptr =
                    game_state.is_lobby ? generate_lobby_message(game_state) : generate_game_message(game_state);
            gui_connection_.write(message_ptr->serialized());
        }
    )
    catch (std::exception &ex) {
        LOG(debug) << "Error in ServerGuiBridge Thread. " << ex.what();
    }

    callback_();
}
