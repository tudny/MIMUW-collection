#include <iostream>
#include "tri_list_concepts.h"
#include "tri_list.h"
//
//int test() {
//
//    std::cout << "Hello, World!" << std::endl;
//
//    auto kw = [](int x) { return x * x; };
//    auto nast = [](int x) { return x + 1; };
//
//    auto kw_nast = compose<int>(kw, nast);
//
//    std::cout << "kw(10)=" << kw(10) << "\n";
//    std::cout << "nast(10)=" << nast(10) << "\n";
//    std::cout << "kw_nast(10)=" << kw_nast(10) << "\n";
//
//    tri_list<int, double, std::string> li;
//
//    std::string s = "Xd";
//
//    std::function<int(int)> f1 = [](int x) -> int { return x + 7; };
//    auto f2 = [](int x) -> int { return x + 5; };
//
//    static_assert(modifier<decltype(f1), int>);
//    static_assert(modifier<decltype(f2), int>);
//
//    auto f3 = compose<int>(f1, f2);
//
//    li.push_back(10);
//    li.push_back(20);
//    li.push_back(20.0);
//    li.push_back(30.0);
//    li.push_back(30);
//    li.push_back(s);
////    li.modify_only<int, decltype([](int x) { return x * x; })>();
//    li.modify_only<int>([](int x) { return x * x + 1; });
//
//    auto li_int = li.range_over<int>();
//
//    for (auto el : li_int) {
//        std::cout << el << "\n";
//    }
//
//
//    tri_list<int, float, bool> l({1, 2, true, 4.f});
//
//    l.modify_only<int>([] (int x) {return x * 2;});
//    l.modify_only<int, decltype([] (int x) {return x + 2;})>();
//
//    for (auto i : l.range_over<float>()) {
//        std::cout << i - std::numbers::pi_v<float> << '\n';
//    }
//
//
//    tri_list<int, double, long long> lii({10, 20, 100l, 12.37});
//
//
//
//    for (auto it = lii.begin(); it != lii.end(); ++it) {
//        std::visit([]<typename T>(T t) {
//            if constexpr (std::is_same_v<T, int>) {
//                std::cout << "int=" << t << "\n";
//            }
//            else if constexpr (std::is_same_v<T, long long>) {
//                std::cout << "ll=" << t << "\n";
//            }
//            else if constexpr (std::is_same_v<T, double>) {
//                std::cout << "double=" << t << "\n";
//            }
//        }, *it);
//    }
//
//
//
//
//    std::vector<int> ints{1, 2, 3, 4};
//
//    auto res = ints | std::views::transform([]<typename T> (T e) { std::cout << e << ":=e\n"; return e * e; });
//    for (auto el : res | std::views::reverse | std::views::filter([]<typename T> (T e) { return e & 1; })) {
//        std::cout << el << "\n";
//    }
//
//    return 0;
//}
//
//int testKomar() {
//    tri_list<int, char, double> l;
//    constexpr int op = int(1e3);
//
//    for (int i = 0; i < op; ++i) {
//        std::cout << i << "\n";
//        l.push_back<int>(0);
//        l.modify_only<int>([](int x) { return x + 1; });
//    }
//
//    for (int i = 0; i < op; ++i) {
//        l.begin();
//        l.end();
//        l.range_over<int>();
//        l.range_over<char>();
//        l.range_over<double>();
//    }
//}
//
//#include <functional>
//
//struct A {};
//
//struct B {};
//
//struct C {};
//
//struct F {
//    A operator()(A x) const { return x; }
//};
//
//C g(C x) { return x; }
//
//int test7() {
//    tri_list<A, B, C> l;
//    l.modify_only<A, F>();
//    l.modify_only<B>([](B x) { return x; });
//    l.modify_only<C>(g);
//}

int main() {
    static_assert(is_tri_list_valid<tri_list<int, float, bool>, int, float, bool>);

    const tri_list<int, bool, double> l{false, 1, 2.0};
    const auto it_beg = l.begin();
    const auto it_end = l.end();
    *it_beg;
    it_beg == it_end;
    it_beg != it_end;

    l.range_over<int>();
    l.range_over<bool>();
    l.range_over<double>();

    *l.range_over<int>().begin();
    *l.range_over<bool>().begin();
    *l.range_over<double>().begin();
}


