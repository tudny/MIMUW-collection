#include "../utils/utils.h"
#include "Connection.h"

#include <utility>

Client::Connection::Connection(
        const std::shared_ptr<tcp::socket> &socket)
        : socket(socket) {}

std::string Client::Connection::read(size_t length) {
    LOG(debug) << "[Reading] from " << get_endpoint();
    size_t read = boost::asio::read(*socket,
                                    boost::asio::buffer(buffer_, length),
                                    boost::asio::transfer_exactly(length));
    assert(read == length);
    return std::string{buffer_.data(), read};
}

void Client::Connection::send(const std::string &message) {
    size_t sent = boost::asio::write(*socket,
                                     boost::asio::buffer(message));

    (void) sent;
    assert(sent == message.size());
}

void Client::Connection::send(const std::shared_ptr<std::string> &message) {
    size_t sent = boost::asio::write(*socket,
                                     boost::asio::buffer(*message));

    (void) sent;
    assert(sent == message->size());
}

void Client::Connection::close() {
    if (socket->is_open()) {
        socket->shutdown(tcp::socket::shutdown_both);
        socket->close();
    }
}
