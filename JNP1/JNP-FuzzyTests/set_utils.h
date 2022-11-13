#include "fuzzy.h"
#include <random>

namespace {
    inline int rd(int l, int r) {
        static std::mt19937 rng(0);
        return std::uniform_int_distribution<int>(l, r)(rng);
    }

    inline TriFuzzyNum random_tri() {
        return TriFuzzyNum(rd(-5, 5), rd(-5, 5), rd(-5, 5));
    }
}

namespace SetUtils {
    inline TriFuzzyNumSet random_set(size_t size) {
        TriFuzzyNumSet set;
        
        for (size_t i = 0; i < size; ++i)
            set.insert(random_tri());
        
        return set;
    }
}
