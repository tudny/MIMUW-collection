#ifndef ROVER_LOCATION_H
#define ROVER_LOCATION_H

#include <cstdint>

using coordinate_t = int64_t;

const static size_t DIRECTIONS_COUNT = 4;

enum Direction {
    NORTH,
    EAST,
    SOUTH,
    WEST
};

inline std::ostream &operator<<(std::ostream &os, const Direction &direction) {
    static std::string names[] = {"NORTH", "EAST", "SOUTH", "WEST"};

    return os << names[direction];
}

struct position_t {
    coordinate_t x, y;

    constexpr position_t operator+(const position_t &other) const {
        return {x + other.x, y + other.y};
    }

    constexpr bool operator!=(const position_t &other) const {
        return (this->x != other.x) || (this->y != other.y);
    }

    inline friend std::ostream &operator<<(std::ostream &os, const position_t &position) {
        return os << "(" << position.x << ", " << position.y << ")";
    }
};

struct location_t {
    position_t position;
    Direction direction;

    friend std::ostream &operator<<(std::ostream &os, const location_t &location) {
        return os << location.position << " " << location.direction;
    }
};

namespace DirectionUtils {

    constexpr inline int get_direction_count() {
        return DIRECTIONS_COUNT;
    }

    inline int value_of_direction(const Direction &direction) {
        return direction;
    }

    inline Direction direction_of_value(const int &value) {
        switch (value) {
            case 0:
                return Direction::NORTH;
            case 1:
                return Direction::EAST;
            case 2:
                return Direction::SOUTH;
            case 3:
                return Direction::WEST;
            default:
                throw std::invalid_argument("No directory of given value");
        }
    }

    inline Direction next_clockwise(const Direction &direction, int delta) {
        int dir_count = get_direction_count();
        int new_value = (((value_of_direction(direction) + delta) % dir_count) + dir_count) % dir_count;
        return direction_of_value(new_value);
    }

    inline Direction next_right(const Direction &direction) {
        return next_clockwise(direction, 1);
    }

    inline Direction next_left(const Direction &direction) {
        return next_clockwise(direction, -1);
    }

    inline Direction opposite_direction(const Direction &direction) {
        return next_clockwise(direction, 2);
    }

    inline position_t get_move_delta(const Direction &direction) {
        static coordinate_t dx[] = {0, 1, 0, -1};
        static coordinate_t dy[] = {1, 0, -1, 0};
        int value = value_of_direction(direction);
        return {dx[value], dy[value]};
    }
}

#endif // ROVER_LOCATION_H
