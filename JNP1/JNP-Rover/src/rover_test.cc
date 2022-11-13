#include <cassert>
#include <sstream>
#include "rover.h"

void assert_equal_output(const auto &rover, std::string_view output) {
    std::stringstream s;
    s << rover;
    assert(s.view() == output);
}


int main() {
    auto builder = RoverBuilder();
    auto forward = move_forward();
    auto rover = RoverBuilder()
            .program_command('A', forward)
            .program_command('B', forward)
            .build();
    rover.land({0, 0}, Direction::NORTH);
    rover.execute("AB");
    assert_equal_output(rover, "(0, 2) NORTH");
}
