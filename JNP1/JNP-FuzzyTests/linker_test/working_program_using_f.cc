#include "../fuzzy.h"
#include "foo_lib_using_fuzzy.h"

#include <iomanip>
#include <iostream>
#include <sstream>

int main() {

    TriFuzzyNum a(2, 4, 8);

    for (real_t i = 0; i <= 10; i += 0.5) {
        size_t spacing = f(a, i) * 20.0;
        std::string spaces(spacing, ' ');
        std::stringstream ss;
        ss << "f(" << i << "): ";
        std::cout << std::setfill(' ') << std::setw(10) << ss.str() << spaces << ".\n";
    }

    return 0;
}
