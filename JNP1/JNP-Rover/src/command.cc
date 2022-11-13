#include "command.h"

location_t MoveForward::evaluate(const location_t &current) {
    auto move_delta = DirectionUtils::get_move_delta(current.direction);
    auto new_position = current.position + move_delta;
    return {new_position, current.direction};
}

location_t MoveBackward::evaluate(const location_t &current) {
    auto move_delta = DirectionUtils::get_move_delta(DirectionUtils::opposite_direction(current.direction));
    auto new_position = current.position + move_delta;
    return {new_position, current.direction};
}

location_t RotateRight::evaluate(const location_t &current) {
    auto new_direction = DirectionUtils::next_right(current.direction);
    return {current.position, new_direction};
}

location_t RotateLeft::evaluate(const location_t &current) {
    auto new_direction = DirectionUtils::next_left(current.direction);
    return {current.position, new_direction};
}

Compose::Compose(std::initializer_list<command_ptr> list) : commands_sequence(list) {}

std::pair<location_t, bool> Compose::execute(const location_t &current) {
    CommandCompositionExecutor executor(current);

    for (auto &command_ptr: commands_sequence) {
        executor.next_command(command_ptr);
        if (executor.is_stopped())
            break;
    }

    return executor.result();
}

location_t Compose::evaluate([[maybe_unused]] const location_t &current) {
    throw std::runtime_error("This command cannot be evaluated.");
}

void Compose::set_sensors(const std::shared_ptr<sensors_collection_t> &collection) {
    Command::set_sensors(collection);

    for (auto &command : commands_sequence) {
        command->set_sensors(collection);
    }
}

command_ptr move_forward() {
    return std::make_shared<MoveForward>();
}

command_ptr move_backward() {
    return std::make_shared<MoveBackward>();
}

command_ptr rotate_right() {
    return std::make_shared<RotateRight>();
}

command_ptr rotate_left() {
    return std::make_shared<RotateLeft>();
}

command_ptr compose(std::initializer_list<command_ptr> list) {
    return std::make_shared<Compose>(list);
}

void CommandCompositionExecutor::next_command(command_ptr &command) {
    if (all_succeeded) {
        auto[new_location, can_move] = command->execute(current_location);
        current_location = new_location;
        all_succeeded &= can_move;
    }
}

std::pair<location_t, bool> CommandCompositionExecutor::result() {
    return {current_location, all_succeeded};
}

void CommandCompositionExecutor::make_stop() {
    all_succeeded = false;
}

bool CommandCompositionExecutor::is_stopped() const {
    return !all_succeeded;
}
