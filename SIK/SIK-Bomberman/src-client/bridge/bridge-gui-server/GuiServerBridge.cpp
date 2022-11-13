/* BOOST */
#include <boost/format.hpp>

/* LOCAL */
#include "GuiServerBridge.h"

GuiServerBridge::GuiServerBridge(
        GUI::Connection &gui_connection,
        Server::Connection &server_connection,
        const callback_t &callback,
        GameState &game_state,
        const std::string &player_name) :
            Bridge(gui_connection,
                   server_connection,
                   callback,
                   game_state),
            player_name(player_name) {}

void GuiServerBridge::operator()() {

    try_forever (
        auto message = gui_connection_.read();
        LOG(trace) << boost::format("[GUI] Received message (%d).") % message.size();

        try {
            auto input_message = InputMessageFactory{}.getInputMessage(message);
            auto client_message_serialized = input_message->generate_client_message()->serialized();

            if (game_state.is_lobby) {
                LOG(debug) << "[GUI] Sending JOIN.";
                JoinClient join_message{player_name};
                server_connection_.write(join_message.serialized());
            }
            else {
                LOG(trace) << boost::format("[GUI] Sending to server (%d).") % client_message_serialized.size();
                server_connection_.write(client_message_serialized);
            }
        }
        catch (DeserializationException &deserializationException) {
            LOG(trace) << "[GUI] Message could not be deserialized.";
        }
    )
    catch (std::exception &ex) {
        LOG(debug) << "Error in GuiServerBridge Thread. " << ex.what();
    }

    callback_();
}
