#include "../fuzzy.h"
#include "../utils.h"
#include "../logger.h"
#include <cassert>

int main() {

    {
        TriFuzzyNum a(10.0, 20.0, 30.0);
        TriFuzzyNum b(20.0, 30.0, 40.0);

        assert((a <=> b) < 0);
        assert(a < b);
    }

    {
        const TriFuzzyNum a(10.0, 20.0, 30.0);
        const TriFuzzyNum b(20.0, 30.0, 40.0);

        assert((a <=> b) < 0);
        assert(a < b);
    }

    {
        TriFuzzyNum a(10.0, 20.0, 30.0);
        TriFuzzyNum b(20.0, 30.0, 40.0);

        assert((b <=> a) > 0);
        assert(b > a);
    }

    {
        const TriFuzzyNum a(10.0, 20.0, 30.0);
        const TriFuzzyNum b(20.0, 30.0, 40.0);

        assert((b <=> a) > 0);
        assert(b > a);
    }

    {
        TriFuzzyNum a(10.0, 20.0, 30.0);
        TriFuzzyNum b(20.0, 30.0, 40.0);

        assert((a <=> b) < 0);
        assert(a < b);
    }

    {
        const TriFuzzyNum a(10.0, 20.0, 30.0);
        const TriFuzzyNum b(20.0, 30.0, 40.0);

        assert((a <=> b) <= 0);
        assert(a <= b);
    }

    {
        TriFuzzyNum a(10.0, 20.0, 30.0);
        TriFuzzyNum b(20.0, 30.0, 40.0);

        assert((b <=> a) >= 0);
        assert(b >= a);
    }

    {
        const TriFuzzyNum a(10.0, 20.0, 30.0);
        const TriFuzzyNum b(20.0, 30.0, 40.0);

        assert((b <=> a) >= 0);
        assert(b >= a);
    }

    {
        TriFuzzyNum a(10.0, 20.0, 30.0);
        TriFuzzyNum b(20.0, 30.0, 40.0);

        assert((a <=> b) != 0);
    }

    {
        const TriFuzzyNum a(10.0, 20.0, 30.0);
        const TriFuzzyNum b(20.0, 30.0, 40.0);

        assert((a <=> b) != 0);
    }

    {
        TriFuzzyNum a(10.0, 20.0, 30.0);
        TriFuzzyNum b(20.0, 30.0, 40.0);

        assert((b <=> a) != 0);
    }

    {
        const TriFuzzyNum a(10.0, 20.0, 30.0);
        const TriFuzzyNum b(20.0, 30.0, 40.0);

        assert((b <=> a) != 0);
    }

    {
        TriFuzzyNum a(10.0, 20.0, 30.0);
        TriFuzzyNum b(20.0, 30.0, 40.0);

        assert(a != b);
    }

    {
        const TriFuzzyNum a(10.0, 20.0, 30.0);
        const TriFuzzyNum b(20.0, 30.0, 40.0);

        assert(a != b);
    }

    {
        TriFuzzyNum a(10.0, 20.0, 30.0);
        TriFuzzyNum b(10.0, 20.0, 30.0);

        assert(a == b);
    }

    {
        const TriFuzzyNum a(10.0, 20.0, 30.0);
        const TriFuzzyNum b(10.0, 20.0, 30.0);

        assert(a == b);
    }

    {
        constexpr TriFuzzyNum a(10.0, 20.0, 30.0);
        constexpr TriFuzzyNum b(20.0, 30.0, 40.0);

        static_assert(a != b);
    }

    {
        constexpr const TriFuzzyNum a(10.0, 20.0, 30.0);
        constexpr const TriFuzzyNum b(20.0, 30.0, 40.0);

        static_assert(a != b);
    }

    {
        constexpr TriFuzzyNum a(10.0, 20.0, 30.0);
        constexpr TriFuzzyNum b(10.0, 20.0, 30.0);

        static_assert(a == b);
    }

    {
        constexpr const TriFuzzyNum a(10.0, 20.0, 30.0);
        constexpr const TriFuzzyNum b(10.0, 20.0, 30.0);

        static_assert(a == b);
    }

    return 0;
}
