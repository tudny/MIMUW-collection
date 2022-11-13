#ifndef SIK_BOMBERMAN_ACCEPTOR_H
#define SIK_BOMBERMAN_ACCEPTOR_H

/* BOOST */
#include <boost/asio.hpp>

/* LOCAL */
#include "../utils/utils.h"

namespace Acceptor {

    using boost::asio::ip::tcp;

    class Acceptor {
    public:
        Acceptor(boost::asio::io_context &ioContext,
                 const uint16_t &port,
                 const callback_t &callback,
                 const new_connection_callback_t &new_connection);

        void operator()();

        void close();

    private:
        boost::asio::io_context &io_context;
        tcp::acceptor acceptor;
        callback_t callback;
        new_connection_callback_t new_connection;
    };

}


#endif //SIK_BOMBERMAN_ACCEPTOR_H
