/* STD */
#include <iostream>

/* BOOST */
#include <boost/thread.hpp>

/* LOCAL */
#include "utils/arguments.h"
#include "gui/Connection.h"
#include "server/Connection.h"
#include "serializer/serializer.h"
#include "bridge/bridge-gui-server/GuiServerBridge.h"
#include "bridge/bridge-server-gui/ServerGuiBridge.h"

int main(int argc, char *argv[]) {
    // Create io_context for ::boost::asio.
    boost::asio::io_context io_context;

    // Load arguments. If error occurs it will call std::exit().
    LOG(debug) << "Parsing arguments";
    auto arguments = Client::parse_arguments(argc, argv);

    try {
        // Create GUI connection.
        LOG(debug) << "Create GUI connection.";
        auto gui_connection = GUI::Connection{
            io_context,
            *arguments->port_,
            *arguments->gui_address_
        };

        // Create Server connection.
        LOG(debug) << "Create Server connection.";
        auto server_connection = Server::Connection{
            io_context,
            *arguments->port_,
            *arguments->server_address_
        };

        // Main thread waking mechanism.
        LOG(debug) << "Main thread waking mechanism.";
        std::mutex mutex;
        std::condition_variable wait_for_bridges;
        std::atomic<std::thread::id> waker;
        auto wake_main_thread = [&]() {
            waker = std::this_thread::get_id();
            wait_for_bridges.notify_all();
        };

        // Creating game state object.
        LOG(debug) << "Creating Game State object.";
        GameState game_state{};

        // Create GUI->SERVER bridge.
        LOG(debug) << "Create GUI->SERVER bridge.";
        auto gui_server_bridge = GuiServerBridge{
                gui_connection,
                server_connection,
                wake_main_thread,
                game_state,
                *arguments->player_name_
        };

        // Create SERVER->GUI bridge.
        LOG(debug) << "Create SERVER->GUI bridge.";
        auto server_gui_bridge = ServerGuiBridge{
                gui_connection,
                server_connection,
                wake_main_thread,
                game_state
        };

        // Run bridges on different threads.
        LOG(debug) << "Run gui_server_bridge_thread.";
        std::thread gui_server_bridge_thread{gui_server_bridge};
        LOG(debug) << "Run server_gui_bridge_thread.";
        std::thread server_gui_bridge_thread{server_gui_bridge};

        // Wait for any thread to finish.
        std::unique_lock<std::mutex> lock{mutex};
        wait_for_bridges.wait(lock, [&] {
            return waker != std::thread::id{};
        });
        LOG(debug) << "Woken up.";

        // Try to shut down sockets gently.
        LOG(debug) << "Trying to close sockets.";
        gui_server_bridge.close();
        LOG(debug) << "Closed gui_server_bridge.";
        server_gui_bridge.close();
        LOG(debug) << "Closed server_gui_bridge.";

        // Wait for threads to end.
        LOG(debug) << "Wait for threads to end.";
        gui_server_bridge_thread.join();
        LOG(debug) << "Joined gui_server_bridge_thread.";
        server_gui_bridge_thread.join();
        LOG(debug) << "Joined server_gui_bridge_thread.";
    }
    catch (std::exception &ex) {
        // If anything fails, end program.
        LOG(debug) << "Some exception was caught in main(). " << ex.what();
        std::exit(EXIT_FAILURE);
    }

    LOG(debug) << "Ending due to exception in thread.";
    std::exit(EXIT_FAILURE);
}
