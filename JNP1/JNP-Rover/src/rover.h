#ifndef ROVER_ROVER_H
#define ROVER_ROVER_H

#include <map>
#include <stdexcept>
#include <memory>

#include "sensor.h"
#include "location.h"
#include "command.h"

using command_shortcut = char;

class RoverNotLanded : public std::runtime_error {
public:
    RoverNotLanded() : std::runtime_error("RoverNotLanded") {};
};

class RoverBuilder;

class Rover {
    friend RoverBuilder;

    using command_mapper_t = std::map<command_shortcut, command_ptr>;

public:
    Rover() : command_mapper(), sensors_collection(std::make_shared<sensors_collection_t>()) {}

    void land(position_t position, Direction direction);

    void execute(const std::string &commands);

    friend std::ostream &operator<<(std::ostream &os, const Rover &rover);

private:
    void add_sensor(sensor_ptr &&sensor);

    void program_command(const command_shortcut &shortcut, command_ptr &&command);
    void program_command(const command_shortcut &shortcut, command_ptr &command);

    command_mapper_t command_mapper;
    std::shared_ptr<sensors_collection_t> sensors_collection;
    bool landed = false, stopped = false;
    location_t location = {{0, 0}, Direction::NORTH};
};

class RoverBuilder {
public:
    RoverBuilder() = default;

    RoverBuilder &program_command(const command_shortcut &shortcut, command_ptr &&command);
    RoverBuilder &program_command(const command_shortcut &shortcut, command_ptr &command);

    RoverBuilder &add_sensor(sensor_ptr &&sensor);

    Rover &&build();

private:
    Rover rover;
};

#endif // ROVER_ROVER_H
