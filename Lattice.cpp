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

static int get_sum_neighbors_general(const LatticeInt &latt, int idx) {
    std::vector<int> coord = latt.indexToCoord(idx);
    std::vector<int> neighbor(coord);
    int sum_neighbors = 0;
    for (int mu = 0; mu < Nd; mu++) {
        neighbor[mu] = cyclic_shift(coord[mu], 1, latt.L);
        sum_neighbors += latt.get_value_const(neighbor);
        neighbor[mu] = cyclic_shift(coord[mu], -1, latt.L);
        sum_neighbors += latt.get_value_const(neighbor);
        neighbor[mu] = coord[mu];
    }

    return sum_neighbors;
}

static int get_sum_neighbors_2d(const LatticeInt &latt, int idx) {
    const int vol = latt.volume();
    const int L = latt.L;
    int up = (idx + L) % vol;
    int down = (idx + vol - latt.L) % vol;
    int left, right;
    if (idx % L == 0) {
        left = idx + L - 1;
    }
    else {
        left = idx - 1;
    }
    if (idx % L == L - 1) {
        right = idx - L + 1;
    }
    else {
        right = idx + 1;
    }

    return latt.values[up] + latt.values[down]
        + latt.values[left] + latt.values[right];
}

int get_sum_neighbors(const LatticeInt &latt, int idx) {
    if (Nd != 2) {
        return get_sum_neighbors_general(latt, idx);
    }
    else {
        return get_sum_neighbors_2d(latt, idx);
    }
}
