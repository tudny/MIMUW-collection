#ifndef SIK_BOMBERMAN_CLIENT_CONNECTION_H
#define SIK_BOMBERMAN_CLIENT_CONNECTION_H


#include <boost/asio.hpp>
#include <boost/array.hpp>
#include "../utils/utils.h"

namespace Client {

    using boost::asio::ip::tcp;

    enum {
        CLIENT_BUFFER_SIZE = 1 << 10
    };

    class Connection {
    public:
        explicit Connection(const std::shared_ptr<tcp::socket> &socket);
        ~Connection() {
            LOG(debug) << "[Connection] destructor called";
        }

        std::string read(size_t length);

        void send(const std::string &message);
        void send(const std::shared_ptr<std::string> &message);

        void close();

        std::string get_endpoint() {
            auto endpoint = socket->remote_endpoint();
            std::stringstream ss;
            ss << endpoint;
            return ss.str();
        }

    private:
        std::shared_ptr<tcp::socket> socket;
        boost::array<char, CLIENT_BUFFER_SIZE> buffer_{};
    };
}


#endif //SIK_BOMBERMAN_CLIENT_CONNECTION_H
