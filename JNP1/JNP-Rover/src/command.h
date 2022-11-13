#ifndef ROVER_COMMAND_H
#define ROVER_COMMAND_H

#include <memory>
#include <set>
#include <vector>
#include "location.h"
#include "sensor.h"

class Command {
    std::shared_ptr<sensors_collection_t> sensors;

    bool check_sensors(const position_t &position) {
        return std::all_of(sensors->begin(), sensors->end(), [&](auto &sensor) {
            return sensor->is_safe(position.x, position.y);
        });
    }

    virtual location_t evaluate(const location_t &current) = 0;

public:

    Command() = default;

    virtual std::pair<location_t, bool> execute(const location_t &current) {
        auto new_location = evaluate(current);
        bool sensors_ok = true;
        if (new_location.position != current.position) {
            sensors_ok = check_sensors(new_location.position);
        }

        return sensors_ok
               ? std::make_pair(new_location, true) : std::make_pair(current, false);
    }

    virtual void set_sensors(const std::shared_ptr<sensors_collection_t> &collection) {
        sensors = collection;
    }
};

using command_ptr = std::shared_ptr<Command>;

class MoveForward : public Command {
    location_t evaluate(const location_t &current) override;
};

command_ptr move_forward();

class MoveBackward : public Command {
    location_t evaluate(const location_t &current) override;
};

command_ptr move_backward();

class RotateRight : public Command {
    location_t evaluate(const location_t &current) override;
};

command_ptr rotate_right();

class RotateLeft : public Command {
    location_t evaluate(const location_t &current) override;
};

command_ptr rotate_left();

class Compose : public Command {
public:
    Compose(std::initializer_list<command_ptr> list);

    std::pair<location_t, bool> execute(const location_t &current) override;

    void set_sensors(const std::shared_ptr<sensors_collection_t> &collection) override;

private:
    location_t evaluate(const location_t &current) override;

private:

    std::vector<command_ptr> commands_sequence;
};

command_ptr compose(std::initializer_list<command_ptr> list);

class CommandCompositionExecutor {
public:
    explicit CommandCompositionExecutor(const location_t &location) : current_location(location) {}

    void next_command(std::shared_ptr<Command> &command);

    std::pair<location_t, bool> result();

    void make_stop();

    bool is_stopped() const;

private:
    location_t current_location;
    bool all_succeeded = true;
};

#endif // ROVER_COMMAND_H
