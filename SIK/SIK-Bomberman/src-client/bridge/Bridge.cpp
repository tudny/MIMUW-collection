/* LOCAL */
#include "Bridge.h"

Bridge::Bridge(
            GUI::Connection &gui_connection,
            Server::Connection &server_connection,
            const callback_t &callback,
            GameState &game_state
        ) :
            gui_connection_(gui_connection),
            server_connection_(server_connection),
            callback_(callback),
            game_state(game_state) {}

void Bridge::close() {
    if (gui_connection_.is_open())
        gui_connection_.close();

    if (server_connection_.is_open())
        server_connection_.close();
}
