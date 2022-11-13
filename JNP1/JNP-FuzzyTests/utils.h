#include "fuzzy.h"

namespace Utils {
    inline bool equals(const TriFuzzyNum &a, const TriFuzzyNum &b) {
        bool equal = true;
        equal &= (a.lower_value() == b.lower_value());
        equal &= (a.modal_value() == b.modal_value());
        equal &= (a.upper_value() == b.upper_value());
        return equal;
    }

    inline bool if_and_only_if(bool p, bool q) {
        return p == q;
    }
}
