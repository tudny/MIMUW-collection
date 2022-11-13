#ifndef SIK_BOMBERMAN_ARGUMENTS_H
#define SIK_BOMBERMAN_ARGUMENTS_H

/* BOOST */
#include <boost/lexical_cast.hpp>

/* STD */
#include <memory>

using port_t = uint16_t;

namespace Server {
    using U8_ptr = std::shared_ptr<uint8_t>;
    using U16_ptr = std::shared_ptr<uint16_t>;
    using U32_ptr = std::shared_ptr<uint32_t>;
    using U64_ptr = std::shared_ptr<uint64_t>;
    using string_ptr = std::shared_ptr<std::string>;

    class Arguments {
    public:

        Arguments(U16_ptr bombTimer,
                  U8_ptr playersCount,
                  U64_ptr turnDuration,
                  U16_ptr explosionRadius,
                  U16_ptr initialBlocks,
                  U16_ptr gameLength,
                  string_ptr serverName,
                  U16_ptr port,
                  U32_ptr seed,
                  U16_ptr sizeX,
                  U16_ptr sizeY);

        const U16_ptr bomb_timer;
        const U8_ptr players_count;
        const U64_ptr turn_duration;
        const U16_ptr explosion_radius;
        const U16_ptr initial_blocks;
        const U16_ptr game_length;
        const string_ptr server_name;
        const U16_ptr port;
        const U32_ptr seed;
        const U16_ptr size_x;
        const U16_ptr size_y;
    };

    class ArgumentsBuilder {
    public:
        ArgumentsBuilder &set_bomb_timer(const uint16_t&);
        ArgumentsBuilder &set_players_count(const uint8_t&);
        ArgumentsBuilder &set_turn_duration(const uint64_t&);
        ArgumentsBuilder &set_explosion_radius(const uint16_t&);
        ArgumentsBuilder &set_initial_blocks(const uint16_t&);
        ArgumentsBuilder &set_game_length(const uint16_t&);
        ArgumentsBuilder &set_server_name(const std::string&);
        ArgumentsBuilder &set_port(const std::uint16_t&);
        ArgumentsBuilder &set_seed(const std::uint32_t&);
        ArgumentsBuilder &set_size_x(const std::uint16_t&);
        ArgumentsBuilder &set_size_y(const std::uint16_t&);


        std::shared_ptr<Arguments> build();

    private:
        U16_ptr bomb_timer;
        U8_ptr players_count;
        U64_ptr turn_duration;
        U16_ptr explosion_radius;
        U16_ptr initial_blocks;
        U16_ptr game_length;
        string_ptr server_name;
        U16_ptr port;
        U32_ptr seed;
        U16_ptr size_x;
        U16_ptr size_y;
    };

    std::shared_ptr<Arguments> parse_arguments(int argc, char *argv[]);
}

#endif //SIK_BOMBERMAN_ARGUMENTS_H
