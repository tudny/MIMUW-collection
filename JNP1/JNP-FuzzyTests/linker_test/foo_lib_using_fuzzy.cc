#include "foo_lib_using_fuzzy.h"
#include "../fuzzy.h"

real_t f(const TriFuzzyNum &num, real_t x) {
    const real_t l = num.lower_value();
    const real_t m = num.modal_value();
    const real_t u = num.upper_value();

    if (l <= x && x <= m) {
        if (l == m) {
            return x == m;
        } else {
            return (x - l) / (m - l);
        }
    } else if (m < x && x <= u) {
        if (u == m) {
            return x == m;
        } else {
            return (u - x) / (u - m);
        }
    } else {
        return 0;
    }
}
