#ifndef SIK_BOMBERMAN_RNG_H
#define SIK_BOMBERMAN_RNG_H

/* STD */
#include <memory>
#include <random>

namespace Server {
    class RNG {
    public:
        explicit RNG(const uint32_t &initial_seed);
        uint32_t value();

    private:
        const uint32_t initial_seed;
        std::minstd_rand rng;
    };

    std::shared_ptr<RNG> make_rng(const uint32_t &seed);
}


#endif //SIK_BOMBERMAN_RNG_H
