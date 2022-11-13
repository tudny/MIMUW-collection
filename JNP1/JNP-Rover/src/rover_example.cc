#include <cassert>
#include <iostream>
#include <memory>
#include "rover.h"

struct TrueSensor : public Sensor {
    bool is_safe([[maybe_unused]] coordinate_t x,
                 [[maybe_unused]] coordinate_t y) override {
        return true;
    }
};

struct FalseSensor : public Sensor {
    bool is_safe([[maybe_unused]] coordinate_t x,
                 [[maybe_unused]] coordinate_t y) override {
        return false;
    }
};

int main() {
    // Podczas budowy łazika można zaprogramawać jego komendy oraz ustawić różne czujniki.
    auto rover = RoverBuilder()
            .program_command('F', move_forward())
            .program_command('B', move_backward())
            .program_command('R', rotate_right())
            .program_command('L', rotate_left())
            .program_command('U', compose({rotate_right(), rotate_right()}))
            .add_sensor(std::make_unique<TrueSensor>())
            .add_sensor(std::make_unique<TrueSensor>())
            .build();

    // Przed lądowaniem nie jest możliwe sterowanie łazikiem.
    std::cout << rover << std::endl; // unknown
    try {
        rover.execute("F");
        assert(false);
    } catch (std::exception const &e) {
    }

    // Po lądowaniu łazik wykonuje przesłane komendy.
    rover.land({0, 0}, Direction::EAST);
    std::cout << rover << std::endl; // (0, 0) EAST
    rover.execute("FFBRLU");
    std::cout << rover << std::endl; // (1, 0) WEST

    // Łazik zatrzymuje się, gdy napotka nieznaną komendę.
    rover.execute("FXFFF");
    std::cout << rover << std::endl; // (0, 0) WEST stopped

    // Łazik wykonuje poprawne komendy.
    rover.execute("FFF");
    std::cout << rover << std::endl; // (-3, 0) WEST

    // Łazik zatrzymuje się również, gdy czujniki zgłaszają niebezpieczeństwo.
    auto broken_rover = RoverBuilder()
            .program_command('X', move_forward())
            .add_sensor(std::make_unique<FalseSensor>())
            .build();
    broken_rover.land({-1, -1}, Direction::WEST);
    broken_rover.execute("X");
    std::cout << broken_rover << std::endl; // (-1, -1) WEST stopped
}
