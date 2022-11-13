#include "../fuzzy.h"
#include "../utils.h"
#include "../logger.h"
#include <cassert>

void _plus() {
    TriFuzzyNum x(10, 20, 30);
    TriFuzzyNum y(40, 50, 60);

    TriFuzzyNum z = x + y;

    assert(Utils::equals(z, TriFuzzyNum(50, 70, 90)));
}

void _pluse() {
    TriFuzzyNum x(10, 20, 30);
    TriFuzzyNum y(40, 50, 60);

    x += y;

    assert(Utils::equals(x, TriFuzzyNum(50, 70, 90)));
}

void _minus() {
    TriFuzzyNum x(10, 20, 30);
    TriFuzzyNum y(40, 50, 60);

    TriFuzzyNum z = x - y;

    assert(Utils::equals(z, TriFuzzyNum(-50, -30, -10)));
}

void _minuse() {
    TriFuzzyNum x(10, 20, 30);
    TriFuzzyNum y(40, 50, 60);

    x -= y;

    assert(Utils::equals(x, TriFuzzyNum(-50, -30, -10)));
}

void _cdot() {
    TriFuzzyNum x(10, 20, 30);
    TriFuzzyNum y(40, 50, 60);

    TriFuzzyNum z = x * y;

    assert(Utils::equals(z, TriFuzzyNum(400, 1000, 1800)));
}

void _cdote() {
    TriFuzzyNum x(10, 20, 30);
    TriFuzzyNum y(40, 50, 60);

    x *= y;

    assert(Utils::equals(x, TriFuzzyNum(400, 1000, 1800)));
}

int main() {

    _plus();
    _pluse();

    _minus();
    _minuse();

    _cdot();
    _cdote();

    return 0;
}
