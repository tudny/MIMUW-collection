#include "RNG.h"

std::shared_ptr<Server::RNG> Server::make_rng(const uint32_t &seed) {
    return std::make_shared<RNG>(seed);
}

Server::RNG::RNG(const uint32_t &initial_seed) : initial_seed(initial_seed), rng(initial_seed) {}

uint32_t Server::RNG::value() {
    return static_cast<uint32_t>(rng());
}
