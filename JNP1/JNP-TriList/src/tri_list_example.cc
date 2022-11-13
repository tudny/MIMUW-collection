#include <vector>
#include <variant>
#include <concepts>
#include <iostream>
#include <numbers>
#include <string>
#include <cassert>
#include <algorithm>
#include <ranges>
#include <limits>
#include <iterator>

#include "tri_list.h"
#include "tri_list_concepts.h"

int main() {
    static_assert(is_tri_list_valid<tri_list<int, float, bool>, int, float, bool>);
    static_assert(modifier<decltype(identity<int>), int>);

    tri_list<int, float, bool> l({1, 2, true, 4.f});

    l.modify_only<int>([] (int x) {return x * 2;});
    l.modify_only<int, decltype([] (int x) {return x + 2;})>();

    std::vector<int> ints;

    for (auto it = l.begin(); it != l.end(); ++it) {
        std::visit([&]<typename T> (T t) -> void {
            if constexpr (std::same_as<T, int>) {
                ints.push_back(t);
            }
            return;
        }, *it);
    }

    assert((std::vector{ints} == std::vector{4, 6}));

    for (auto i : l.range_over<float>()) {
        std::cout << i - std::numbers::pi_v<float> << '\n';
    }

    using namespace std::string_literals;
    tri_list<int, std::string, std::vector<std::string>> l2 = {
            1, 2, "This", 21, "is",
            std::vector{"This"s, "is"s, "a"s, "vector"s},
            "a", 3, "string", 720
    };

    std::ranges::copy(l2.range_over<std::string>(), std::ostream_iterator<std::string>(std::cout, ", "));
    std::cout << '\n';

    auto bad_iter = std::ranges::max_element(l.range_over<int>());
    static_assert(std::same_as<decltype(bad_iter), std::ranges::dangling>);

    using V = std::variant<int, std::string, std::vector<std::string>>;

    auto good_iter = std::ranges::max_element(l2, std::ranges::less{}, [](const V& v) -> int {
        if (std::holds_alternative<int>(v)) {
            return std::get<int>(v);
        } else {
            return std::numeric_limits<int>::min();
        }
    });

    assert((std::get<int>(*good_iter) == 720));

    tri_list<int, float, bool> l3 = {true, 4.f, true, 5.f, true, 1, 1, 1};

    l3.modify_only<bool>([]([[maybe_unused]] bool b) {
        return false;
    });

    std::cout << std::ranges::count(l3.range_over<int>(), 1) << '\n';

    float acc = 0;
    for (const auto& var: l3) {
        if (std::holds_alternative<float>(var)) {
            acc += std::get<float>(var);
        }
        if (std::holds_alternative<int>(var)) {
            acc += std::get<int>(var);
        }
    }

    std::cout << acc << '\n';

    tri_list<char, bool, int> l4 = {true};

    for (int i = 0; i < 2000; ++i) {
        auto f = compose<int>(
                [](int c) {return c + 1;},
                identity<int>);

        l4.modify_only<int>(f);
    }

    l4.push_back(0);

    assert(std::get<int>(*++std::begin(l4)) == 2000);

    l4.reset<int>();

    auto view_copy = l4.range_over<int>();

    assert(*std::ranges::begin(view_copy) == 0);
}
