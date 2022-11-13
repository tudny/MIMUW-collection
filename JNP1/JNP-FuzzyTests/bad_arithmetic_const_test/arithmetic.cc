#include "../fuzzy.h"
#include "../utils.h"
#include "../logger.h"
#include <cassert>

void plus() {
    const TriFuzzyNum x(10, 20, 30);
    const TriFuzzyNum y(40, 50, 60);

    TriFuzzyNum z = x + y;

    assert(Utils::equals(z, TriFuzzyNum(50, 70, 90)));
}

void pluse() {
    const TriFuzzyNum x(10, 20, 30);
    const TriFuzzyNum y(40, 50, 60);

    x += y;

    assert(Utils::equals(x, TriFuzzyNum(50, 70, 90)));
}

void minus() {
    const TriFuzzyNum x(10, 20, 30);
    const TriFuzzyNum y(40, 50, 60);

    TriFuzzyNum z = x - y;

    assert(Utils::equals(z, TriFuzzyNum(-50, -30, -10)));
}

void minuse() {
    const TriFuzzyNum x(10, 20, 30);
    const TriFuzzyNum y(40, 50, 60);

    x -= y;

    assert(Utils::equals(x, TriFuzzyNum(-50, -30, -10)));
}

void cdot() {
    const TriFuzzyNum x(10, 20, 30);
    const TriFuzzyNum y(40, 50, 60);

    TriFuzzyNum z = x * y;

    assert(Utils::equals(z, TriFuzzyNum(400, 1000, 1800)));
}

void cdote() {
    const TriFuzzyNum x(10, 20, 30);
    const TriFuzzyNum y(40, 50, 60);

    x *= y;

    assert(Utils::equals(x, TriFuzzyNum(400, 1000, 1800)));
}

int main() {

    plus();
    pluse();

    minus();
    minuse();

    cdot();
    cdote();

    return 0;
}
