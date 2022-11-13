#include "Sender.h"

#include <utility>

Client::Sender::Sender(std::shared_ptr<Client::Connection> connection, std::shared_ptr<PlayerInstance> instance)
        : connection(std::move(connection)), instance(std::move(instance)) {}

void Client::Sender::operator()() {
    std::mutex mutex;
    std::condition_variable wait_for_messages;
    std::atomic_bool waken{false};
    instance->output_actions.set_listener_on_change([&]() {
        waken = true;
        wait_for_messages.notify_all();
    });

    try_forever (
        // Is there are no messages to send we sleep for 1 second. Otherwise, we send them.
        if (instance->output_actions.empty()) {
            std::unique_lock lock{mutex};
            wait_for_messages.wait_for(lock, std::chrono::seconds(1), [&]() {
                // This predicate is being checked before falling asleep.
                return waken == true;
            });
        }

        while (!instance->output_actions.empty()) {
            LOG(debug) << "Sender. Something to send.";
            auto message = instance->output_actions.take();
            if (message) {
                connection->send(message->serialize());
                LOG(debug) << "Sending.";
            }
        }
    )
    catch (...) {
        LOG(debug) << "Error in Sender.";
    }

    instance->disconnected = true;
}
