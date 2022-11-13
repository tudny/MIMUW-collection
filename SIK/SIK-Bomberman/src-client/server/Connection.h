#ifndef SIK_BOMBERMAN_SERVER_CONNECTION_H
#define SIK_BOMBERMAN_SERVER_CONNECTION_H

/* BOOST */
#include <boost/asio.hpp>
#include <boost/array.hpp>

/* LOCAL */
#include "../utils/arguments.h"

namespace Server {

    namespace asio = boost::asio;
    using boost::asio::ip::tcp;

    const size_t SERVER_BUFFER_SIZE = 1 << 18; ///< Size of a buffer.

    /**
     * This class represent the connection to the Server.
     */
    class Connection {
    public:
        /**
         * The constructor of a Connection. It opens socket.
         * @param io_context the io context
         * @param port the port to send from
         * @param gui_address the address of a Server
         */
        Connection(
                asio::io_context &io_context,
                const port_t &port,
                Client::Address &address);

        /**
         * Closes the socket.
         */
        void close();

        /**
         * Reads the given amount of bytes.
         * @param length the amount of bytes to read
         * @return the message represented as an array of bytes (std::string).
         */
        std::string read(size_t length);

        /**
         * Sends datagram to Server.
         * @param message the message to be sent.
         */
        void write(const std::string &message);

        /**
         * Checks if the socket is open.
         * @return is socket open
         */
        bool is_open();

    private:

        /**
         * Tries connect to all endpoints with given port.
         * @param port the port
         */
        void try_all_endpoints_from(const port_t &port);

        /**
         * Sets all socket options.
         * NO_DELAY
         * REUSE_ADDRESS
         * is protocol of an endpoint is ipv6 sets v6_only=false
         * @param endpoint endpoint with protocol type
         */
        void set_socket_options(const tcp::endpoint &endpoint);

        asio::ip::tcp::resolver::iterator server_endpoints_;
        tcp::socket socket_;
        boost::array<char, SERVER_BUFFER_SIZE> buffer_{};
    };
}


#endif //SIK_BOMBERMAN_SERVER_CONNECTION_H
