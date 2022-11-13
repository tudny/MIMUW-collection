#ifndef ROBOTS_ARGUMENTS_H
#define ROBOTS_ARGUMENTS_H

/* BOOST */
#include <boost/lexical_cast.hpp>

using port_t = uint16_t;

namespace Client {
    struct Address {
        std::string address;
        port_t port;

        explicit Address(const std::string &str);

        friend std::ostream &operator<<(std::ostream &os, const Client::Address &_address) {
            return os << _address.address << ":" << _address.port;
        }
    };

    using port_ptr = std::shared_ptr<port_t>;
    using string_ptr = std::shared_ptr<std::string>;
    using address_ptr = std::shared_ptr<Address>;

    class Arguments {
    public:
        Arguments(address_ptr guiAddress,
                  string_ptr playerName,
                  port_ptr port,
                  address_ptr serverAddress);

    public:
        const address_ptr gui_address_;
        const string_ptr player_name_;
        const port_ptr port_;
        const address_ptr server_address_;
    };

    class ArgumentsBuilder {
    public:
        ArgumentsBuilder &set_gui_address(const std::string&);
        ArgumentsBuilder &set_player_name(const std::string&);
        ArgumentsBuilder &set_port(const port_t&);
        ArgumentsBuilder &set_server_address(const std::string&);

        std::shared_ptr<Arguments> build();
    private:
        address_ptr _gui_address;
        string_ptr _player_name;
        port_ptr _port;
        address_ptr _server_address;
    };

    std::shared_ptr<Arguments> parse_arguments(int argc, char *argv[]);
}

#endif //ROBOTS_ARGUMENTS_H
