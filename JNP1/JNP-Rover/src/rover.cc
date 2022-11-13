#include "rover.h"

RoverBuilder &RoverBuilder::program_command(const command_shortcut &shortcut, command_ptr &&command) {
    rover.program_command(shortcut, std::move(command));
    return *this;
}

RoverBuilder &RoverBuilder::program_command(const command_shortcut &shortcut, command_ptr &command) {
    rover.program_command(shortcut, command);
    return *this;
}

RoverBuilder &RoverBuilder::add_sensor(sensor_ptr &&sensor) {
    rover.add_sensor(std::move(sensor));
    return *this;
}

Rover &&RoverBuilder::build() {
    return std::move(rover);
}

void Rover::add_sensor(sensor_ptr &&sensor) {
    sensors_collection->insert(std::move(sensor));
}

void Rover::program_command(const command_shortcut &shortcut, command_ptr &&command) {
    auto command_it = command_mapper.find(shortcut);
    if (command_it != command_mapper.end()) {
        command_mapper.erase(command_it);
    }

    command->set_sensors(sensors_collection);
    command_mapper.insert(std::make_pair(shortcut, command));
}

void Rover::program_command(const command_shortcut &shortcut, command_ptr &command) {
    command->set_sensors(sensors_collection);
    command_mapper[shortcut] = command;
}

void Rover::land(position_t position, Direction direction) {
    location = {position, direction};
    landed = true;
    stopped = false;
}

void Rover::execute(const std::string &commands) {
    if (!landed)
        throw RoverNotLanded();

    stopped = false;

    CommandCompositionExecutor executor(location);

    for (const auto &shortcut: commands) {
        auto command_it = command_mapper.find(shortcut);

        if (command_it != command_mapper.end())
            executor.next_command(command_it->second);
        else
            executor.make_stop();

        if (executor.is_stopped())
            break;
    }

    auto[new_location, all_success] = executor.result();
    location = new_location;
    if (!all_success)
        stopped = true;
}

std::ostream &operator<<(std::ostream &os, const Rover &rover) {
    if (!rover.landed)
        return os << "unknown";

    os << rover.location;

    if (rover.stopped)
        os << " stopped";

    return os;
}
