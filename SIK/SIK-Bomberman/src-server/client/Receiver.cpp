#include "Receiver.h"

#include <utility>

Client::Receiver::Receiver(const std::shared_ptr<Client::Connection> &connection, const std::shared_ptr<PlayerInstance> &instance)
        : connection(connection), instance(instance) {}

void Client::Receiver::operator()() {
    try_forever (
        auto client_message = Models::ClientMessageFactory(*connection);
        instance->input_actions.push(client_message);
        LOG(trace) << "Received";
    )
    catch (std::exception &ex) {
        LOG(error) << "error receiver " << ex.what();
    }

    instance->disconnected = true;
}
