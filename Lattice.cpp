#include "Lattice.h"

// Returns a random choice between +1 and -1
static int random_spin() {
    thread_local static std::default_random_engine generator;
    thread_local static std::bernoulli_distribution distribution(0.5);
    return distribution(generator) * 2 - 1;
}

LatticeInt random_lattice(int L) {
    LatticeInt ret(L);
    for (int idx = 0; idx < ret.volume(); idx++) {
        ret.values[idx] = random_spin();
    }
    return ret;
}
