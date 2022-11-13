/* STD */
#include <string>
#include <memory>
#include <iostream>
#include <utility>

/* LOCAL */
#include "arguments.h"

/* BOOST */
#include <boost/program_options.hpp>

Client::Arguments::Arguments(address_ptr guiAddress,
                             string_ptr playerName,
                             port_ptr port,
                             address_ptr serverAddress) :
                                 gui_address_(std::move(guiAddress)),
                                 player_name_(std::move(playerName)),
                                 port_(std::move(port)),
                                 server_address_(std::move(serverAddress)) {}

Client::ArgumentsBuilder &Client::ArgumentsBuilder::set_gui_address(const std::string &gui_address) {
    this->_gui_address = std::make_shared<Address>(gui_address);
    return *this;
}

Client::ArgumentsBuilder &Client::ArgumentsBuilder::set_player_name(const std::string &player_name) {
    this->_player_name = std::make_shared<std::string>(player_name);
    return *this;
}

Client::ArgumentsBuilder &Client::ArgumentsBuilder::set_port(const port_t &port) {
    this->_port = std::make_shared<port_t>(port);
    return *this;
}

Client::ArgumentsBuilder &Client::ArgumentsBuilder::set_server_address(const std::string &server_address) {
    this->_server_address = std::make_shared<Address>(server_address);
    return *this;
}

std::shared_ptr<Client::Arguments> Client::ArgumentsBuilder::build() {
    return std::make_shared<Client::Arguments>(
            Client::Arguments(_gui_address,
                              _player_name,
                              _port,
                              _server_address));
}

namespace po = boost::program_options;

static void print_usage(const po::options_description &description, int exit_code) {
    std::cout << description << std::endl;
    std::exit(exit_code);
}

static void print_error(const std::string &message) {
    std::cerr << message << std::endl;
}

#define GUI_ADDRESS "gui-address"
#define HELP "help"
#define PLAYER_NAME "player-name"
#define PORT "port"
#define SERVER_ADDRESS "server-address"
#define ADDITIONAL "additional"

std::shared_ptr<Client::Arguments> Client::parse_arguments(int argc, char *argv[]) {
    po::options_description description("Bomberman client.\nUsage: "
                                            + std::string(argv[0]) + " <options>\nThese are options");
    po::positional_options_description positional_options_description;
    positional_options_description.add(ADDITIONAL, -1);
    description.add_options()
            (GUI_ADDRESS",d", po::value<std::string>()->required(),
             "<(nazwa hosta):(port) lub (IPv4):(port) lub (IPv6):(port)>")
            (HELP",h",
                    "Print help information")
            (PLAYER_NAME",n", po::value<std::string>()->required(),
                    "<String>")
            (PORT",p", po::value<port_t>()->required())
            (SERVER_ADDRESS",s", po::value<std::string>()->required(),
                    "<(nazwa hosta):(port) lub (IPv4):(port) lub (IPv6):(port)>");

    po::variables_map variables_map;

    try {
        po::store(po::command_line_parser(argc, argv)
                        .options(description)
                        .positional(positional_options_description)
                        .run(),
                    variables_map);

        if (variables_map.count(HELP)) {
            print_usage(description, 0);
            __builtin_unreachable();
        }

        if (variables_map.contains(ADDITIONAL)) {
            throw std::runtime_error("What is '" + variables_map[ADDITIONAL].as<std::string>() + "'?");
        }

        po::notify(variables_map);

        return ArgumentsBuilder{}
                .set_gui_address(variables_map[GUI_ADDRESS].as<std::string>())
                .set_player_name(variables_map[PLAYER_NAME].as<std::string>())
                .set_port(variables_map[PORT].as<port_t>())
                .set_server_address(variables_map[SERVER_ADDRESS].as<std::string>())
                .build();
    }
    catch (const std::exception &exception) {
        print_error("Wrong usage! " + std::string(exception.what()));
        print_usage(description, 1);
    }
    __builtin_unreachable();
}

Client::Address::Address(const std::string &str) {
    auto div = str.find_last_of(':');
    address = str.substr(0, div);
    port = boost::lexical_cast<port_t>(str.substr(div + 1));
}
