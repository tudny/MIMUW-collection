#include "../fuzzy.h"
#include "../logger.h"
#include "../utils.h"
#include <cassert>

int main() {

    constexpr TriFuzzyNum num(10.0, 20.0, 30.0);

    constexpr real_t l_c = num.lower_value();
    constexpr real_t m_c = num.modal_value();
    constexpr real_t u_c = num.upper_value();

    static_assert(l_c == 10.0);
    static_assert(m_c == 20.0);
    static_assert(u_c == 30.0);

    real_t l = num.lower_value();
    real_t m = num.modal_value();
    real_t u = num.upper_value();

    assert(l == l_c);
    assert(m == m_c);
    assert(u == u_c);

    Logger::log("Constexpr check.");

    return 0;
}
