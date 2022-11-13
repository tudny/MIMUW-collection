/* STD */
#include <chrono>
#include <iostream>
#include <utility>

/* LOCAL */
#include "arguments.h"

/* BOOST */
#include <boost/program_options.hpp>

Server::Arguments::Arguments(Server::U16_ptr bombTimer,
                             Server::U8_ptr playersCount,
                             Server::U64_ptr turnDuration,
                             Server::U16_ptr explosionRadius,
                             Server::U16_ptr initialBlocks,
                             Server::U16_ptr gameLength,
                             Server::string_ptr serverName,
                             Server::U16_ptr port,
                             Server::U32_ptr seed,
                             Server::U16_ptr sizeX,
                             Server::U16_ptr sizeY)
        : bomb_timer(std::move(bombTimer)),
            players_count(std::move(playersCount)),
            turn_duration(std::move(turnDuration)),
            explosion_radius(std::move(explosionRadius)),
            initial_blocks(std::move(initialBlocks)),
            game_length(std::move(gameLength)),
            server_name(std::move(serverName)),
            port(std::move(port)), seed(std::move(seed)),
            size_x(std::move(sizeX)), size_y(std::move(sizeY)) {}

Server::ArgumentsBuilder &Server::ArgumentsBuilder::set_bomb_timer(const uint16_t &bomb_timer_) {
    bomb_timer = std::make_shared<uint16_t>(bomb_timer_);
    return *this;
}

Server::ArgumentsBuilder &Server::ArgumentsBuilder::set_players_count(const uint8_t &players_count_) {
    players_count = std::make_shared<uint8_t>(players_count_);
    return *this;
}

Server::ArgumentsBuilder &Server::ArgumentsBuilder::set_turn_duration(const uint64_t &turn_duration_) {
    turn_duration = std::make_shared<uint64_t>(turn_duration_);
    return *this;
}

Server::ArgumentsBuilder &Server::ArgumentsBuilder::set_explosion_radius(const uint16_t &explosion_radius_) {
    explosion_radius = std::make_shared<uint16_t>(explosion_radius_);
    return *this;
}

Server::ArgumentsBuilder &Server::ArgumentsBuilder::set_initial_blocks(const uint16_t &initial_blocks_) {
    initial_blocks = std::make_shared<uint16_t>(initial_blocks_);
    return *this;
}

Server::ArgumentsBuilder &Server::ArgumentsBuilder::set_game_length(const uint16_t &game_length_) {
    game_length = std::make_shared<uint16_t>(game_length_);
    return *this;
}

Server::ArgumentsBuilder &Server::ArgumentsBuilder::set_server_name(const std::string &server_name_) {
    server_name = std::make_shared<std::string>(server_name_);
    return *this;
}

Server::ArgumentsBuilder &Server::ArgumentsBuilder::set_port(const uint16_t &port_) {
    port = std::make_shared<uint16_t>(port_);
    return *this;
}

Server::ArgumentsBuilder &Server::ArgumentsBuilder::set_seed(const uint32_t &seed_) {
    seed = std::make_shared<uint32_t>(seed_);
    return *this;
}

Server::ArgumentsBuilder &Server::ArgumentsBuilder::set_size_x(const uint16_t &size_x_) {
    size_x = std::make_shared<uint16_t>(size_x_);
    return *this;
}

Server::ArgumentsBuilder &Server::ArgumentsBuilder::set_size_y(const uint16_t &size_y_) {
    size_y = std::make_shared<uint16_t>(size_y_);
    return *this;
}

static uint32_t default_seed() {
    return static_cast<uint32_t>(
            std::chrono::system_clock::now().time_since_epoch().count());
}

std::shared_ptr<Server::Arguments> Server::ArgumentsBuilder::build() {
    return std::make_shared<Server::Arguments>(
            bomb_timer,
            players_count,
            turn_duration,
            explosion_radius,
            initial_blocks,
            game_length,
            server_name,
            port,
            seed,
            size_x,
            size_y);
}

#define BOMB_TIMER "bomb-timer"
#define PLAYERS_COUNT "players-count"
#define TURN_DURATION "turn-duration"
#define EXPLOSION_RADIUS "explosion-radius"
#define HELP "help"
#define INITIAL_BLOCKS "initial-blocks"
#define GAME_LENGTH "game-length"
#define SERVER_NAME "server-name"
#define PORT "port"
#define SEED "seed"
#define SIZE_X "size-x"
#define SIZE_Y "size-y"
#define ADDITIONAL "additional"

namespace po = boost::program_options;

static void print_usage(const po::options_description &description, int exit_code) {
    std::cout << description << std::endl;
    std::exit(exit_code);
}

static void print_error(const std::string &message) {
    std::cerr << message << std::endl;
}

std::shared_ptr<Server::Arguments> Server::parse_arguments(int argc, char **argv) {
    po::options_description description("Bomberman server.\nUsage: "
                            + std::string(argv[0]) + " <options>\nThese are options");
    po::positional_options_description positional_options_description;
    positional_options_description.add(ADDITIONAL, -1);
    description.add_options()
            (BOMB_TIMER",b", po::value<uint16_t>()->required())
            (PLAYERS_COUNT",c", po::value<uint16_t>()->required())
            (TURN_DURATION",d", po::value<uint64_t>()->required(), "milisecunds")
            (EXPLOSION_RADIUS",e", po::value<uint16_t>()->required())
            (HELP",h", "Print help information")
            (INITIAL_BLOCKS",k", po::value<uint16_t>()->required())
            (GAME_LENGTH",l", po::value<uint16_t>()->required())
            (SERVER_NAME",n", po::value<std::string>()->required())
            (PORT",p", po::value<uint16_t>()->required())
            (SEED",s", po::value<uint32_t>()->default_value(default_seed()), "optional")
            (SIZE_X",x", po::value<uint16_t>()->required())
            (SIZE_Y",y", po::value<uint16_t>()->required());

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
        auto player_count = variables_map[PLAYERS_COUNT].as<uint16_t>();
        if (player_count > UINT8_MAX)
            throw std::runtime_error("Players count is not a UINT8.");

        return ArgumentsBuilder{}
                .set_bomb_timer(variables_map[BOMB_TIMER].as<uint16_t>())
                .set_players_count(static_cast<uint8_t>(player_count))
                .set_turn_duration(variables_map[TURN_DURATION].as<uint64_t>())
                .set_explosion_radius(variables_map[EXPLOSION_RADIUS].as<uint16_t>())
                .set_initial_blocks(variables_map[INITIAL_BLOCKS].as<uint16_t>())
                .set_game_length(variables_map[GAME_LENGTH].as<uint16_t>())
                .set_server_name(variables_map[SERVER_NAME].as<std::string>())
                .set_port(variables_map[PORT].as<uint16_t>())
                .set_seed(variables_map[SEED].as<uint32_t>())
                .set_size_x(variables_map[SIZE_X].as<uint16_t>())
                .set_size_y(variables_map[SIZE_Y].as<uint16_t>())
                .build();
    }
    catch (const std::exception &exception) {
        print_error("Wrong usage! " + std::string(exception.what()));
        print_usage(description, 1);
    }
    __builtin_unreachable();
}
