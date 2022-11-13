#include "../fuzzy.h"
#include "../utils.h"
#include "../logger.h"
#include <sstream>
#include <cassert>

int main() {

    std::stringstream ss;

    TriFuzzyNum num1(1.25, 2.25, 3.25);
    TriFuzzyNum num2(2.25, 3.25, 1.25);
    constexpr TriFuzzyNum num3(1.25, 2.25, 3.25);
    constinit static TriFuzzyNum num4(1.25, 2.25, 3.25);

    ss << num1 << " " << num2 << " " << num3 << " " << num4;

    std::string result = ss.str();
    std::string expected = "(1.25, 2.25, 3.25) (1.25, 2.25, 3.25) (1.25, 2.25, 3.25) (1.25, 2.25, 3.25)";
    assert(result.compare(expected) == 0);

    return 0;
}
