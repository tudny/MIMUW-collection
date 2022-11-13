#include <vector>
#include <concepts>
#include <memory>
#include <functional>

#include "tri_list_concepts.h"

template<typename T>
T identity(const T &t) {
    return t;
}

template<typename T, modifier<T> F1, modifier<T> F2>
modifier<T> auto compose(F1 f1, F2 f2) {
    return [=](const T &t) mutable {
        return std::invoke(f1, std::invoke(f2, t));
    };
}

template<typename... T>
constexpr bool exactly_one(T... predicates) {
    return (... + predicates) == 1;
}

template<typename T1, typename T2, typename T3, typename T>
concept proper_type = exactly_one(std::same_as<T, T1>, std::same_as<T, T2>, std::same_as<T, T3>);

template<typename T1, typename T2, typename T3>
class tri_list {
    using type_t = std::variant<T1, T2, T3>;
    using container_t = std::vector<type_t>;

public:
    tri_list() = default;

    tri_list(std::initializer_list<type_t> list) : container_original(list) {}

    tri_list(const tri_list &other) = default;

    tri_list(tri_list &&other) = default;

    template<typename T>
    requires proper_type<T1, T2, T3, T>
    void push_back(const T &t) {
        container_original.push_back(t);
    }

    template<typename T, modifier<T> F>
    requires proper_type<T1, T2, T3, T>
    void modify_only(F m = F{}) {
        if constexpr (std::is_same_v<T, T1>) {
            t1_modifiers = compose<T1>(m, t1_modifiers);
        }
        else if constexpr (std::is_same_v<T, T2>) {
            t2_modifiers = compose<T2>(m, t2_modifiers);
        }
        else if constexpr (std::is_same_v<T, T3>) {
            t3_modifiers = compose<T3>(m, t3_modifiers);
        }
    }

    template<typename T>
    requires proper_type<T1, T2, T3, T>
    void reset() {
        if constexpr (std::is_same_v<T, T1>) {
            t1_modifiers = identity<T1>;
        }
        else if constexpr (std::is_same_v<T, T2>) {
            t2_modifiers = identity<T2>;
        }
        else if constexpr (std::is_same_v<T, T3>) {
            t3_modifiers = identity<T3>;
        }
    }

    template<typename T>
    requires proper_type<T1, T2, T3, T>
    auto range_over() const {
        auto result = container_original | std::views::filter([](type_t ele) {
            return std::holds_alternative<T>(ele);
        }) | std::views::transform([](type_t ele) {
            return std::get<T>(ele);
        });

        if constexpr (std::is_same_v<T, T1>) {
            return result | std::views::transform(t1_modifiers);
        }
        else if constexpr (std::is_same_v<T, T2>) {
            return result | std::views::transform(t2_modifiers);
        }
        else if constexpr (std::is_same_v<T, T3>) {
            return result | std::views::transform(t3_modifiers);
        }
    }

    class iterator {
        typename container_t::const_iterator it;

        std::function<T1(T1)> t1_modifiers = identity<T1>;
        std::function<T2(T2)> t2_modifiers = identity<T2>;
        std::function<T3(T3)> t3_modifiers = identity<T3>;

    public:
        using iterator_category = std::bidirectional_iterator_tag;
        using difference_type = std::iter_difference_t<typename container_t::const_iterator>;
        using value_type = type_t;
        using pointer = type_t *;
        using reference = type_t;

        iterator() noexcept = default;

        explicit iterator(typename container_t::const_iterator _it,
                          std::function<T1(T1)> t1_modifiers = identity<T1>,
                          std::function<T2(T2)> t2_modifiers = identity<T2>,
                          std::function<T3(T3)> t3_modifiers = identity<T3>) :
                it(_it),
                t1_modifiers(t1_modifiers),
                t2_modifiers(t2_modifiers),
                t3_modifiers(t3_modifiers) {}

        iterator operator++(int) noexcept {
            iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        iterator operator--(int) noexcept {
            iterator tmp = *this;
            --(*this);
            return tmp;
        }

        iterator &operator++() noexcept {
            it++;
            return *this;
        }

        iterator &operator--() noexcept {
            it--;
            return *this;
        }

        reference operator*() const {
            auto ele = *it;

            if (std::holds_alternative<T1>(ele)) {
                return type_t(t1_modifiers(std::get<T1>(ele)));
            }
            else if (std::holds_alternative<T2>(ele)) {
                return type_t(t2_modifiers(std::get<T2>(ele)));
            }
            else if (std::holds_alternative<T3>(ele)) {
                return type_t(t3_modifiers(std::get<T3>(ele)));
            }

            throw std::runtime_error("No such type");
        }

        pointer *operator->() const {
            return std::make_shared<type_t>(this->operator*());
        }

        bool operator==(const iterator &other) const {
            return it == other.it;
        }

        bool operator!=(const iterator &other) const {
            return it != other.it;
        }
    };

    auto begin() const {
        return iterator(container_original.begin(), t1_modifiers, t2_modifiers, t3_modifiers);
    }

    auto end() const {
        return iterator(container_original.end(), t1_modifiers, t2_modifiers, t3_modifiers);;
    }

private:
    container_t container_original;
    std::function<T1(T1)> t1_modifiers = identity<T1>;
    std::function<T2(T2)> t2_modifiers = identity<T2>;
    std::function<T3(T3)> t3_modifiers = identity<T3>;
};
