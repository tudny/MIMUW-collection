/* LOCAL */
#include "Connection.h"
#include "../utils/utils.h"

namespace GUI {

    /**
     * Get address in format from arguments and create ::asio::udp::endpoints.
     * @param gui_address address to parse
     * @param io_context the io_context
     * @return endpoints basen on address
     */
    udp::resolver::iterator parse_endpoint(const Client::Address &gui_address, asio::io_context &io_context) {
        udp::resolver resolver{io_context};
        udp::resolver::query query{gui_address.address, std::to_string(gui_address.port)};
        return resolver.resolve(query);
    }

    Connection::Connection(
            asio::io_context &io_context,
                const port_t &port,
                const Client::Address &gui_address
            ) :
               socket_{
                 io_context,
                 udp::endpoint{udp::v6(), port}
               },
               gui_endpoints_(parse_endpoint(gui_address, io_context)) {}

    void Connection::close() {
        LOG(trace) << "Closing gui (UDP) socket.";
        boost::system::error_code error_code;
        socket_.shutdown(boost::asio::socket_base::shutdown_both, error_code);
        socket_.close(error_code);
        LOG(trace) << "Close(UDP)@ErrorCode=" << error_code;
        LOG(trace) << "Closed gui (UDP) socket.";
    }

    std::string Connection::read() {
        udp::endpoint remote_endpoint;
        size_t length = socket_.receive_from(asio::buffer(recv_buffer), remote_endpoint);

        // https://github.com/agluszak/mimuw-sik-2022-public
        // P: Czy możemy być pewni, że wiadomość od GUI przyszła z podanego adresu i wiadomości do GUI są wysyłane
        //    z podanego portu? Innymi słowy, czy wiadomości od GUI mamy odbierać przez receive, czy receive_from
        //    (i analogicznie wysyłać przez send, czy send_to)?
        // O: Adres i port GUI, które podaje się w kliencie, służą do wysyłania wiadomości od klienta do GUI.
        //    GUI może wysyłać komunikaty z portów efemerycznych. Ale ogólnie najlepiej nic nie zakładać o adresie
        //    GUI i być gotowym na odbieranie (poprawnych) wiadomości od kogokolwiek
        // assert(gui_endpoint == remote_endpoint);

        return std::string{recv_buffer.data(), length};
    }

    bool Connection::try_sending(const udp::endpoint &endpoint, const std::string &message) {
        try {
            size_t len = socket_.send_to(asio::buffer(message), endpoint);
            if (len != message.size())
                return false;
        }
        catch (...) {
            return false;
        }

        return true;
    }

    void Connection::write(const std::string &message) {
        size_t successful_sends = 0;
        udp::resolver::iterator end;
        udp::resolver::iterator begin = gui_endpoints_;

        // We try to send datagram to any endpoint resolved by the resolver.
        // If all failed we crash the program.
        while (begin != end) {
            if (successful_sends += try_sending(*begin++, message))
                break;
        }

        if (!successful_sends)
            throw NoEndpointException{};
    }

    bool Connection::is_open() {
        return socket_.is_open();
    }
}

