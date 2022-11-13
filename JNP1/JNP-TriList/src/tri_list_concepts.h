#ifndef TRI_LIST_CONCEPTS_H
#define TRI_LIST_CONCEPTS_H

#include <concepts>
#include <ranges>
#include <variant>

template <typename T1, typename T2, typename T3>
class tri_list;

// implementation provided to students
template <typename F, typename T>
concept modifier = requires (F f, T t) {
    requires std::invocable<F, T>;
    {f(t)} -> std::convertible_to<T>;
};

// implementation provided to students
template <typename TriList, typename T1, typename T2, typename T3>
concept is_tri_list_valid = requires (TriList l) {
    {l} -> std::convertible_to<tri_list<T1, T2, T3>>;
    {l} -> std::ranges::viewable_range;
    {*l.begin()} -> std::same_as<std::variant<T1, T2, T3>>;
    {l.template range_over<T1>()} -> std::ranges::view;
};

#endif
