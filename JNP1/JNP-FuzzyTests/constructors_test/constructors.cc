#include "../fuzzy.h"
#include "../logger.h"
#include "../utils.h"
#include <cassert>

TriFuzzyNum num1(1.25, 2.25, 3.25);
TriFuzzyNum num2(2.25, 3.25, 1.25);
constexpr TriFuzzyNum num3(1.25, 2.25, 3.25);
constinit static TriFuzzyNum num4(1.25, 2.25, 3.25);

void constructors() {
    using std::cout;

    cout << num1 << "\n";
    cout << num2 << "\n";
    cout << num3 << "\n";
    cout << num4 << "\n";

    assert(Utils::equals(num1, num2));
    assert(Utils::equals(num1, num3));
    assert(Utils::equals(num1, num4));
    assert(Utils::equals(num2, num3));
    assert(Utils::equals(num2, num4));
    assert(Utils::equals(num3, num4));

    Logger::log("TriFuzzyNum real_t constructor seems to work.");
}

void coping_constructor() {
    const TriFuzzyNum a(1, 2, 3);
    TriFuzzyNum b(a);

    assert(Utils::equals(a, b));

    TriFuzzyNum x(4, 5, 6);
    TriFuzzyNum y(x);

    assert(Utils::equals(x, y));

    constexpr TriFuzzyNum n(4, 5, 6);
    TriFuzzyNum m(n);

    assert(Utils::equals(n, m));

    constexpr static TriFuzzyNum l(4, 5, 6);
    TriFuzzyNum k(l);

    assert(Utils::equals(l, k));

    Logger::log("TriFuzzyNum coping constructor seems to work.");
}

void moving_constructor() {
    TriFuzzyNum a(TriFuzzyNum(10.0, 20.0, 30.0));

    assert(Utils::equals(a, TriFuzzyNum(10.0, 20.0, 30.0)));

    TriFuzzyNum b(30, 40, 50);
    TriFuzzyNum c(std::move(b));

    assert(Utils::equals(c, TriFuzzyNum(30, 40, 50)));

    Logger::log("TriFuzzyNum moving constructor seems to work.");
}

int main() {
    constructors();
    coping_constructor();
    moving_constructor();

    return 0;
}
