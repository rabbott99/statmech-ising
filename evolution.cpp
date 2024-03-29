#include "Lattice.h"
#include "evolution.h"

#include <cmath>
#include <random>

static double uniform_random() {
    thread_local static std::default_random_engine generator;
    thread_local static std::uniform_real_distribution<double> distribution(0,1);
    return distribution(generator);
}

void run_sweep(LatticeInt &latt, double beta, double muB) {
    assert(beta > 0);
    for (int idx = 0; idx < latt.volume(); idx++) {
        int sum_neighbors = get_sum_neighbors(latt, idx);

        // If s_i = +1, then \delta s_i = -2
        // If s_i = -1, then \delta s_i = +2
        int delta_s_i = -2 * latt.values[idx];
        double deltaH = -delta_s_i * sum_neighbors - muB * delta_s_i;

        // Accept-reject step
        if (deltaH > 0) {
            // Note that e^{-\beta \Delta H} < 1 if and only if \Delta H > 0,
            // so we only need to run a check if deltaH > 0
            double threshold = exp(- beta * deltaH);
            double roll = uniform_random();
            if (roll > threshold) {
                // Reject
                continue;
            }
        }
        // Accept
        latt.values[idx] = -latt.values[idx];
    }
}
