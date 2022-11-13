#ifndef FUZZY_H
#define FUZZY_H

#include <iostream>
#include <set>
#include <array>
#include <cmath>
#include <vector>
#include <tuple>
#include <numeric>
#include <unordered_set>

using std::move;
using std::forward;
using std::tie;
using std::array;
using std::tuple_cat;
using std::sort;

using real_t = double;

class TriFuzzyNum {

    static const size_t FUZZY_NUMBER_ELEMENTS = 3;
    // std::array can be compared with standard <=> operator lexicographically
    using rang_t = array<real_t, FUZZY_NUMBER_ELEMENTS>;

    real_t l, m, u;

    constexpr rang_t rang() const noexcept;
    constexpr void normalize_order() noexcept;
    static constexpr void sort_three_values(real_t&, real_t&, real_t&) noexcept;

    // Method for performing operation *f : (real_t, real_t) -> real_t on every pair of corresponding elements.
    TriFuzzyNum &math_operator(const TriFuzzyNum &other, void (*f)(real_t&, real_t));

    TriFuzzyNum operator- () const;

    // Returns a tuple of form { l, m, u } for operator==
    constexpr auto to_tuple() const;

public:
    constexpr TriFuzzyNum(const real_t &x, const real_t &y, const real_t &z) noexcept;

    // Copy constructor
    TriFuzzyNum(const TriFuzzyNum &other) noexcept;

    // Move constructor
    TriFuzzyNum(TriFuzzyNum &&other) noexcept;

    // Getters
    constexpr real_t lower_value() const;
    constexpr real_t modal_value() const;
    constexpr real_t upper_value() const;

    // Copy assignment
    TriFuzzyNum &operator= (const TriFuzzyNum&) = default;

    // Move assignment
    TriFuzzyNum &operator= (TriFuzzyNum&&) = default;

    // Simple maths operators
    TriFuzzyNum operator+ (const TriFuzzyNum &other) const;
    TriFuzzyNum operator- (const TriFuzzyNum &other) const;
    TriFuzzyNum operator* (const TriFuzzyNum &other) const;

    // Changing maths operators
    TriFuzzyNum &operator+= (const TriFuzzyNum &other);
    TriFuzzyNum &operator-= (const TriFuzzyNum &other);
    TriFuzzyNum &operator*= (const TriFuzzyNum &other);

    // Comparison operators
    constexpr bool operator== (const TriFuzzyNum &rhs) const;
    constexpr bool operator!= (const TriFuzzyNum &rhs) const;

    std::partial_ordering operator<=> (const TriFuzzyNum &rhs) const;
    bool operator<  (const TriFuzzyNum &rhs) const;
    bool operator<= (const TriFuzzyNum &rhs) const;
    bool operator>  (const TriFuzzyNum &rhs) const;
    bool operator>= (const TriFuzzyNum &rhs) const;
};

/* ************************************ */
/* TriFuzzyNum constexpr implementation */
/* ************************************ */

constexpr TriFuzzyNum::rang_t TriFuzzyNum::rang() const noexcept {
    double z = (u - l) + sqrt(1 + (u - m) * (u - m)) + sqrt(1 + (m - l) * (m - l));
    double x = ((u - l) * m + sqrt(1 + (u - m) * (u - m)) * l + sqrt(1 + (m - l) * (m - l)) * u) / z;
    double y = (u - l) / z;

    return { x - y / 2, 1 - y, m };
}

// Sorting values a, b, c in minimal number of comparisons.
constexpr void TriFuzzyNum::sort_three_values(real_t &a, real_t &b, real_t &c) noexcept {
    if (a > c) std::swap(a, c);
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
}

constexpr void TriFuzzyNum::normalize_order() noexcept {
    sort_three_values(this->l, this->m, this->u);
}

constexpr TriFuzzyNum::TriFuzzyNum(const real_t &x, const real_t &y, const real_t &z) noexcept : l(x), m(y), u(z) {
    normalize_order();
}

constexpr real_t TriFuzzyNum::lower_value() const {
    return l;
}

constexpr real_t TriFuzzyNum::modal_value() const {
    return m;
}

constexpr real_t TriFuzzyNum::upper_value() const {
    return u;
}

constexpr auto TriFuzzyNum::to_tuple() const {
    return std::make_tuple(l, m, u);
}

constexpr bool TriFuzzyNum::operator==(const TriFuzzyNum &rhs) const {
    return this->to_tuple() == rhs.to_tuple();
}

constexpr bool TriFuzzyNum::operator!=(const TriFuzzyNum &rhs) const {
    return !(*this == rhs);
}

struct TriFuzzyNumHasher {
    size_t operator()(const TriFuzzyNum &) const;
};

std::ostream &operator<<(std::ostream &os, const TriFuzzyNum &triFuzzyNum);

/* **************************************** */
/* TriFuzzyNum constexpr implementation end */
/* **************************************** */

class TriFuzzyNumSet {
    using tri_fuzzy_num_set_t = std::unordered_multiset<TriFuzzyNum, TriFuzzyNumHasher>;

    real_t l_sum = 0, m_sum = 0, u_sum = 0;
    tri_fuzzy_num_set_t tri_fuzzy_set;

    void init();

    void count_single(const TriFuzzyNum &num);

    void decount_single(const TriFuzzyNum &num);

    void impact_single(const TriFuzzyNum &num, real_t (*f)(real_t, real_t));

public:

    // Constructors
    TriFuzzyNumSet() = default;

    // Coping constructor
    TriFuzzyNumSet(const TriFuzzyNumSet &other) = default;

    // Moving constructor
    TriFuzzyNumSet(TriFuzzyNumSet &&other) = default;

    // List initializer constructor
    TriFuzzyNumSet(std::initializer_list<TriFuzzyNum> list);

    // Copy assignment
    TriFuzzyNumSet &operator= (const TriFuzzyNumSet&) = default;

    // Move assignment
    TriFuzzyNumSet &operator= (TriFuzzyNumSet&&) = default;

    // Insertion
    void insert(const TriFuzzyNum &element);
    void insert(TriFuzzyNum &&element);

    // Remove
    void remove(const TriFuzzyNum &element);

    // Arithmetic mean
    TriFuzzyNum arithmetic_mean() const;
};

/* ***************** */
/* TriFuzzyNum utils */
/* ***************** */

consteval TriFuzzyNum crisp_number(real_t v) noexcept {
    return TriFuzzyNum(v, v, v);
}

inline constinit const TriFuzzyNum crisp_zero = crisp_number(0);

/* ********************* */
/* TriFuzzyNum utils end */
/* ********************* */


#endif //FUZZY_H
