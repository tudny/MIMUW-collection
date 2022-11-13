#ifndef SIK_BOMBERMAN_GUI_CONNECTION_H
#define SIK_BOMBERMAN_GUI_CONNECTION_H

/* BOOST */
#include <boost/asio.hpp>

/* LOCAL */
#include "../utils/arguments.h"

namespace GUI {

    namespace asio = boost::asio;
    using boost::asio::ip::udp;

    /**
     * An exception thrown when no datagram was sent.
     */
    class NoEndpointException : public std::exception {};

    const size_t GUI_BUFFER_SIZE = 1 << 4; ///< Size of a buffer. No more is needed

    /**
     * This class represent the connection to the GUI.
     */
    class Connection {
    public:
        /**
         * The constructor of a Connection. It opens socket.
         * @param io_context the io context
         * @param port the port to send from
         * @param gui_address the address of a GUI
         */
        Connection(
                asio::io_context &io_context,
                const port_t &port,
                const Client::Address &gui_address);

        /**
         * Closes the socket.
         */
        void close();

        /**
         * Reads the datagram.
         * @return the datagram represented as an array of bytes (std::string).
         */
        std::string read();

        /**
         * Sends datagram to GUI.
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
         * Tries sending the datagram to an endpoint.
         * @param endpoint the endpoint
         * @param message the datagram
         * @return was the sent successful
         */
        bool try_sending(const udp::endpoint &endpoint, const std::string &message);

        udp::socket socket_; ///< connection socket
        udp::resolver::iterator gui_endpoints_; ///< GUI endpoint
        boost::array<char, GUI_BUFFER_SIZE> recv_buffer{}; ///< buffer for messages
    };
}


#endif //SIK_BOMBERMAN_GUI_CONNECTION_H
