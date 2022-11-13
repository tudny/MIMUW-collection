#include "../fuzzy.h"
#include "../utils.h"
#include "../logger.h"
#include <cassert>

consteval bool foo() {
    TriFuzzyNum zero(0.0, 0.0, 0.0);
    TriFuzzyNum not_zero(1.0, 0.0, 0.0);
    return zero.upper_value() != not_zero.upper_value() && !(zero == not_zero);
}

const constinit bool foo_var = foo();

int main() {

    std::cout << foo_var << "\n";
    std::cout << foo() << "\n";

    assert(foo());
    assert(foo_var);
    static_assert(foo());
    static_assert(foo_var);

    return 0;
}
