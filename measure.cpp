#include "measure.h"
#include "evolution.h"

// Computes -\sum_{<ij>} s_i s_j
static double spin_energy(const LatticeInt &latt) {
    int ret = 0;
    for (int idx = 0; idx < latt.volume(); idx++) {
        std::vector<int> coord = latt.indexToCoord(idx);
        std::vector<int> neighbor(coord);
        int sum_neighbors = get_sum_neighbors(latt, idx);
        ret += latt.values[idx] * sum_neighbors;
    }
    // The factor of 1/2 comes from double counting
    return -ret / 2;
}

double magnetization(const LatticeInt &latt) {
    int total = 0;
    for (int idx = 0; idx < latt.volume(); idx++) {
        total += latt.values[idx];
    }
    return (double) total / latt.volume();
}

MeasurementResults measure(const LatticeInt &latt, double muB) {
    MeasurementResults ret;
    ret.magnetization = magnetization(latt);
    const double M = ret.magnetization;
    ret.magnetization_squared = M * M;

    // H = -J \sum_<ij> s_i s_j - mu B \sum_i s_i with J = 1
    ret.energy = spin_energy(latt) - muB * ret.magnetization;
    ret.energy_squared = ret.energy * ret.energy;

    // By choosing a sign for the magnetization at 0 field we prevent the
    // positive and negative magnetization states from cancelling out and
    // producing 0
    ret.magnetization = std::abs(ret.magnetization);
    return ret;
}

std::vector<MeasurementResults> measurement_run(LatticeInt &latt,
        double beta, double muB,
        int thermalization_trajectories, int measurement_separation,
        int num_measurements) {

    auto run_sweeps = [&](int n) {
        for (int i = 0; i < n; i++) {
            run_sweep(latt, beta, muB);
        }
    };

    run_sweeps(thermalization_trajectories);

    std::vector<MeasurementResults> ret;
    ret.reserve(num_measurements);
    for (int i = 0; i < num_measurements; i++) {
        run_sweeps(measurement_separation);
        ret.push_back(measure(latt, muB));
    }
    return ret;
}

// Assumes the measurements have been averaged over several trajectories
PhysicalResults physical_results(const MeasurementResults &measurements,
        int lattice_volume, double T) {
    PhysicalResults ret;
    ret.energy = measurements.energy;

    const double H = measurements.energy;
    const double H2 = measurements.energy_squared;
    ret.heat_capacity = (H2 - H * H) / (lattice_volume * T * T);

    const double M = measurements.magnetization;
    const double M2 = measurements.magnetization_squared;
    ret.magnetization = M;
    ret.susceptibility = (M2 - M * M) * lattice_volume / T;

    return ret;
}
