#include <iostream>
#include <array>
#include <cassert>
#include <cstdlib>

constexpr const size_t N = 100'001;

template <typename Iterator, typename K>
constexpr void fill_constexpr(Iterator begin, Iterator end, K key) {
    while (begin != end) {
        *(begin++) = key;
    }
}

constexpr auto calculate_sieve() {
    std::array<bool, N> sieve{false};
    fill_constexpr(sieve.begin(), sieve.end(), true);

    sieve[0] = sieve[1] = false;

    for (size_t p = 2; p * p < N; ++p) {
        if (sieve[p]) {
            for (size_t i = p * p; i < N; i += p) {
                sieve[i] = false;
            }
        }
    }
    
    return sieve;
}

constexpr const auto is_prime = calculate_sieve();

int main() {

    int a;
    std::cin >> a;

    std::cout << "Is " << a << " prime?\n" 
              << std::boolalpha << is_prime[a] << "\n";

    return 0;
}

// int main(int argc, char *argv[]) {

//     assert(argc > 1);
//     int a = atoi(argv[1]);

//     std::cout << "Is " << a << " prime?\n" 
//               << std::boolalpha << is_prime(a) << "\n";

//     return 0;
// }
