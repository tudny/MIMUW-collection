#include "../fuzzy.h"
#include "../logger.h"
#include "../utils.h"
#include <cassert>

void coping_assignment() {
    const TriFuzzyNum a(1, 2, 3);
    TriFuzzyNum b(4, 5, 6);
    b = a;

    assert(Utils::equals(a, b));

    TriFuzzyNum x(4, 5, 6);
    TriFuzzyNum y(7, 9, 9);
    y = x;

    assert(Utils::equals(x, y));

    constexpr TriFuzzyNum n(4, 5, 6);
    TriFuzzyNum m(9, 10, 11);
    m = n;

    assert(Utils::equals(n, m));

    constexpr static TriFuzzyNum l(4, 5, 6);
    TriFuzzyNum k(12, 13, 14);
    k = l;

    assert(Utils::equals(l, k));

    Logger::log("TriFuzzyNum coping assignment seems to work.");
}

void moving_assignment() {
    TriFuzzyNum a(30, 40, 50);
    a = TriFuzzyNum(10.0, 20.0, 30.0);

    assert(Utils::equals(a, TriFuzzyNum(10.0, 20.0, 30.0)));

    TriFuzzyNum b(30, 40, 50);
    TriFuzzyNum c(70, 80, 90);
    c = std::move(b);

    assert(Utils::equals(c, TriFuzzyNum(30, 40, 50)));

    Logger::log("TriFuzzyNum moving assignment seems to work.");
}

int main() {
    coping_assignment();
    moving_assignment();

    return 0;
}
