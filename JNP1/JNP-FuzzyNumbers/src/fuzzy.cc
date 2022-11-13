#include "fuzzy.h"
#include <functional>

/*  ***********  */
/*  TriFuzzyNum  */
/*  ***********  */

TriFuzzyNum::TriFuzzyNum(const TriFuzzyNum &other) noexcept : l(other.l), m(other.m), u(other.u) {}

TriFuzzyNum::TriFuzzyNum(TriFuzzyNum &&other) noexcept : l(move(other.l)), m(move(other.m)), u(move(other.u)) {}

TriFuzzyNum &TriFuzzyNum::math_operator(const TriFuzzyNum &other, void (*f)(real_t&, real_t)) {
    (*f)(this->l, other.l);
    (*f)(this->m, other.m);
    (*f)(this->u, other.u);

    normalize_order();

    return *this;
}

TriFuzzyNum TriFuzzyNum::operator+(const TriFuzzyNum &other) const {
    return TriFuzzyNum(*this) += other;
}

TriFuzzyNum TriFuzzyNum::operator-(const TriFuzzyNum &other) const {
    return TriFuzzyNum(*this) -= other;
}

TriFuzzyNum TriFuzzyNum::operator*(const TriFuzzyNum &other) const {
    return TriFuzzyNum(*this) *= other;
}

TriFuzzyNum TriFuzzyNum::operator-() const {
    return TriFuzzyNum(-this->l, -this->m, -this->u);
}

TriFuzzyNum &TriFuzzyNum::operator+=(const TriFuzzyNum &other) {
    return math_operator(other, [](real_t &x, real_t y) { x += y; });
}

TriFuzzyNum &TriFuzzyNum::operator-=(const TriFuzzyNum &other) {
    return math_operator(-other, [](real_t &x, real_t y) { x += y; });
}

TriFuzzyNum &TriFuzzyNum::operator*=(const TriFuzzyNum &other) {
    return math_operator(other, [](real_t &x, real_t y) { x *= y; });
}

std::ostream &operator<<(std::ostream &os, const TriFuzzyNum &triFuzzyNum) {
    return os << "("
              << triFuzzyNum.lower_value() << ", "
              << triFuzzyNum.modal_value() << ", "
              << triFuzzyNum.upper_value() << ")";
}

std::partial_ordering TriFuzzyNum::operator<=>(const TriFuzzyNum &rhs) const {
    auto lhs_rang = this->rang();
    auto rhs_rang = rhs.rang();

    return lhs_rang <=> rhs_rang;
}

bool TriFuzzyNum::operator<(const TriFuzzyNum &rhs) const {
    return (*this <=> rhs) < 0;
}

bool TriFuzzyNum::operator<=(const TriFuzzyNum &rhs) const {
    return *this < rhs || *this == rhs;
}

bool TriFuzzyNum::operator>(const TriFuzzyNum &rhs) const {
    return (*this <=> rhs) > 0;
}

bool TriFuzzyNum::operator>=(const TriFuzzyNum &rhs) const {
    return *this > rhs || *this == rhs;
}

size_t TriFuzzyNumHasher::operator()(const TriFuzzyNum &num) const {
    using std::hash;

    return (hash<real_t>()(num.lower_value()) << 1) ^
           ((hash<real_t>()(num.modal_value()) << 1) >> 1) ^
           hash<real_t>()(num.upper_value());
}

/*  **************  */
/*  TriFuzzyNumSet  */
/*  **************  */

void TriFuzzyNumSet::init() {
    for (const TriFuzzyNum &ele : tri_fuzzy_set) {
        count_single(ele);
    }
}

void TriFuzzyNumSet::count_single(const TriFuzzyNum &num) {
    impact_single(num, [](real_t x, real_t y) { return x + y; });
}

void TriFuzzyNumSet::decount_single(const TriFuzzyNum &num) {
    impact_single(num, [](real_t x, real_t y) { return x - y; });
}

void TriFuzzyNumSet::impact_single(const TriFuzzyNum &num, real_t (*f)(real_t, real_t)) {
    l_sum = (*f)(l_sum, num.lower_value());
    m_sum = (*f)(m_sum, num.modal_value());
    u_sum = (*f)(u_sum, num.upper_value());
}

TriFuzzyNumSet::TriFuzzyNumSet(std::initializer_list<TriFuzzyNum> list) : tri_fuzzy_set(list) {
    init();
}

void TriFuzzyNumSet::insert(const TriFuzzyNum &element) {
    count_single(element);
    tri_fuzzy_set.insert(element);
}

void TriFuzzyNumSet::insert(TriFuzzyNum &&element) {
    count_single(element);
    tri_fuzzy_set.insert(move(element));
}

void TriFuzzyNumSet::remove(const TriFuzzyNum &element) {
    auto element_iterator = tri_fuzzy_set.find(element);
    if (element_iterator != tri_fuzzy_set.end()) {
        decount_single(element);
        tri_fuzzy_set.erase(element_iterator);
    }
}

TriFuzzyNum TriFuzzyNumSet::arithmetic_mean() const {
    if (tri_fuzzy_set.empty())
        throw std::length_error("TriFuzzyNumSet::arithmetic_mean - the set is empty.");

    auto div = static_cast<double>(tri_fuzzy_set.size());
    return TriFuzzyNum(l_sum / div, m_sum / div, u_sum / div);
}


