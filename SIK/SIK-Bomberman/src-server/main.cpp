/* STD */
#include <iostream>

/* LOCAL */
#include "utils/arguments.h"
#include "utils/RNG.h"
#include "utils/utils.h"
#include "client/Acceptor.h"
#include "threads/set.h"
#include "client/Connection.h"
#include "client/Sender.h"
#include "client/Receiver.h"
#include "game/GameLoop.h"
#include "client/PlayerInstance.h"

/* BOOST */
#include <boost/asio.hpp>
#include <boost/endian/conversion.hpp>
#include <thread>

int main(int argc, char *argv[]) {
    // Create io_context for ::boost::asio.
    boost::asio::io_context io_context;

    // Load arguments. If error occurs it will call std::exit().
    LOG(debug) << "Parsing arguments";
    auto arguments = Server::parse_arguments(argc, argv);
    LOG(debug) << "Initializing RNG";
    auto rng = Server::make_rng(*arguments->seed);

    // Main thread waking mechanism.
    LOG(debug) << "Main thread waking mechanism.";
    std::mutex mutex;
    std::condition_variable wait_for_bridges;
    std::atomic<std::thread::id> waker;
    std::atomic_bool waken_up = false;
    auto wake_main_thread = [&] {
        waken_up = true;
        waker = std::this_thread::get_id();
        wait_for_bridges.notify_all();
    };

    // Collection of connections.
    safe::set<std::shared_ptr<Client::Connection>> connections;
    std::shared_ptr<safe::set<std::shared_ptr<PlayerInstance>>> instances =
            std::make_shared<safe::set<std::shared_ptr<PlayerInstance>>>();
    std::set<std::shared_ptr<std::thread>> workers;
    safe::queue<std::shared_ptr<Models::ServerMessage>> current_messages;
    std::map<std::shared_ptr<PlayerInstance>, std::shared_ptr<Client::Sender>> senders;
    std::map<std::shared_ptr<PlayerInstance>, std::shared_ptr<Client::Receiver>> receivers;

    auto hello = std::make_shared<Models::Hello>(
            Models::String{*arguments->server_name},
            Models::UINT8{*arguments->players_count},
            Models::UINT16{*arguments->size_x},
            Models::UINT16{*arguments->size_y},
            Models::UINT16{*arguments->game_length},
            Models::UINT16{*arguments->explosion_radius},
            Models::UINT16{*arguments->bomb_timer}
    );

    auto new_connection = [&](const std::shared_ptr<boost::asio::ip::tcp::socket> &socket) {
        if (waken_up) return;
        LOG(debug) << "New connection.";
        auto connection = std::make_shared<Client::Connection>(socket);
        connections.add(connection);
        auto player_instance = std::make_shared<PlayerInstance>();
        instances->add(player_instance);
        player_instance->input_actions.set_listener_on_change([&] { instances->on_change(); });
        auto sender = std::make_shared<Client::Sender>(connection, player_instance);
        senders[player_instance] = sender;
        auto receiver = std::make_shared<Client::Receiver>(connection, player_instance);
        receivers[player_instance] = receiver;
        std::thread(std::ref(*sender)).detach();
        std::thread(std::ref(*receiver)).detach();
        player_instance->output_actions.push(hello);
        current_messages.perform([&](std::queue<std::shared_ptr<Models::ServerMessage>> &queue) {
            auto q_copy = queue;
            while (!q_copy.empty()) {
                player_instance->output_actions.push(q_copy.front());
                q_copy.pop();
            }
        });
    };

    LOG(debug) << "Creating game loop.";
    auto game_loop = Game::GameLoop{ arguments, instances, rng, current_messages, senders, receivers, connections };
    LOG(debug) << "Start game loop thread.";
    std::shared_ptr<std::thread> game_loop_thread = std::make_shared<std::thread>(std::ref(game_loop));
    workers.insert(game_loop_thread);

    LOG(debug) << "Creating acceptor.";
    auto acceptor = Acceptor::Acceptor{
            io_context,
            *arguments->port,
            wake_main_thread,
            new_connection
    };

    LOG(debug) << "Starting acceptor_thread.";
    std::shared_ptr<std::thread> acceptor_thread = std::make_shared<std::thread>(std::ref(acceptor));
    workers.insert(acceptor_thread);

    // Wait for any thread to finish.
    std::unique_lock<std::mutex> lock{mutex};
    wait_for_bridges.wait(lock, [&] {
        return waker != std::thread::id{};
    });
    waken_up = true;
    LOG(debug) << "Woken up.";

    LOG(debug) << "Interrupting threads.";
    ignore_exception ( acceptor.close(); )
    connections.for_each([](auto &connection) {
        ignore_exception (
            connection->close();
        )
    });

    // Wait for threads to end.
    LOG(debug) << "Wait for threads to end.";
    for (auto &worker : workers) {
        LOG(trace) << "Joining thread=0x'" << std::hex << std::hash<std::thread::id>{}(worker->get_id()) << "'.";
        ignore_exception ( worker->join(); )
    }
    LOG(debug) << "Joined acceptor_thread.";

    return 0;
}
