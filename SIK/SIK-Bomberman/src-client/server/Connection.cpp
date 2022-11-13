/* LOCAL */
#include "Connection.h"
#include "../utils/utils.h"

/* STD */
#include <iostream>

/* BOOST */
#include <boost/format.hpp>

namespace Server {

    static asio::ip::tcp::resolver::iterator parse_endpoints(asio::io_context &io_context, Client::Address &address) {
        asio::ip::tcp::resolver resolver(io_context);
        asio::ip::tcp::resolver::query query{address.address, std::to_string(address.port)};
        asio::ip::tcp::resolver::iterator iter = resolver.resolve(query);

        return iter;
    }

    Connection::Connection(
            asio::io_context &io_context,
            const port_t &port,
            Client::Address &address) :
               server_endpoints_{
                 parse_endpoints(
                   io_context,
                   address
                 )
               },
               socket_{io_context} {

        try_all_endpoints_from(port);
    }

    void Connection::set_socket_options(const tcp::endpoint &endpoint) {
        socket_.set_option(asio::ip::tcp::no_delay{true});
        if (endpoint.protocol() == tcp::v6())
            socket_.set_option(asio::ip::v6_only{false});
        socket_.set_option(asio::socket_base::reuse_address{true});
    }

    void Connection::close() {
        socket_.close();
    }

    std::string Connection::read(size_t length) {
        size_t read = asio::read(
                socket_,
                asio::buffer(buffer_, length),
                boost::asio::transfer_exactly(length));

        if  (read != length)
            throw std::runtime_error{"Error reading from the socket."};

        return std::string{buffer_.data(), read};
    }

    void Connection::write(const std::string &message) {
        size_t sent = asio::write(
                socket_,
                asio::buffer(message));

        if (message.size() != sent)
            throw std::runtime_error{"Error sending via socket."};
    }

    bool Connection::is_open() {
        return socket_.is_open();
    }

    void Connection::try_all_endpoints_from(const port_t &port) {
        asio::ip::tcp::resolver::iterator end;
        boost::system::error_code ec;

        LOG(debug) << "Iterative endpoints.";

        while (server_endpoints_ != end) {
            tcp::endpoint server_endpoint_ = *server_endpoints_;
            LOG(debug) << boost::format("Trying new endpoint: '%s'") % server_endpoint_;

            if (socket_.is_open())
                socket_.close(ec);

            socket_.open(server_endpoint_.protocol(), ec);
            set_socket_options(server_endpoint_);
            socket_.bind(
                    tcp::endpoint{
                            server_endpoint_.protocol(),
                            port
                    });
            socket_.connect(server_endpoint_, ec);
            if (!ec) {
                LOG(debug) << (boost::format("Successfully connected to: '%s'") % server_endpoint_);
                break;
            }

            ++server_endpoints_;
        }

        if (ec) {
            LOG(debug) << boost::format("An error occurred during connecting.");
            throw boost::system::system_error(ec);
        }
    }
}
