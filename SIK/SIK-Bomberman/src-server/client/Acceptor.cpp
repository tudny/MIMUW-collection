/* LOCAL */
#include "Acceptor.h"

/* STD */
#include <iostream>

Acceptor::Acceptor::Acceptor(boost::asio::io_context &ioContext,
                             const uint16_t &port,
                             const callback_t &callback,
                             const new_connection_callback_t &new_connection)
        : io_context(ioContext),
          acceptor(io_context, tcp::endpoint(tcp::v6(), port)),
          callback(callback),
          new_connection(new_connection) {}

void Acceptor::Acceptor::operator()() {
    try_forever(
        std::shared_ptr<tcp::socket> socket = std::make_shared<tcp::socket>(io_context);
        acceptor.accept(*socket);
        socket->set_option(tcp::no_delay{true});
        new_connection(socket);
    )
    catch (std::exception &ex) {
        LOG(debug) << "Acceptor failed. " << ex.what();
    }

    callback();
}

void Acceptor::Acceptor::close() {
    if (acceptor.is_open()) {
        acceptor.close();
    }
}
