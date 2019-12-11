#include <iostream>

#include "Lattice.h"
#include "evolution.h"
#include "measure.h"

template <typename T>
std::vector<T> jackknife_resample(const std::vector<T> &samples) {
    assert(samples.size() > 1);
    T total = samples[0];
    for (int i = 1; i < samples.size(); i++) total += samples[i];

    std::vector<T> ret(samples.size(), total);
    for (int i = 0; i < samples.size(); i++) {
        ret[i] -= samples[i];
        ret[i] *= 1.0 / (samples.size() - 1);
    }

    return ret;
}

void print_result(std::ostream &os, const PhysicalResults &res) {
    os << res.energy << " " << res.heat_capacity << " "
        << res.susceptibility;
}

void print_jackknife(std::ostream &os, const std::vector<PhysicalResults> &results) {
    for (int i = 0; i < results.size(); i++) {
        print_result(os, results[i]);
        if (i < results.size() - 1) {
            os << " ";
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: ./ising L T" << std::endl;
        exit(1);
    }
    const int L = atoi(argv[1]);
    double T = atof(argv[2]);
    double beta = 1.0 / T;

    const int thermalization_iterations = 1e5;
    const int measurement_separation = 10;
    const int num_measurements = 3e4;

    LatticeInt start = random_lattice(L);
    std::vector<MeasurementResults> measurements = measurement_run(start,
            beta, 0.0, thermalization_iterations,
            measurement_separation, num_measurements);

    std::vector<MeasurementResults> jackknife = jackknife_resample(measurements);
    std::vector<PhysicalResults> results;
    results.reserve(jackknife.size());
    for (auto sample: jackknife) {
        results.push_back(physical_results(sample, start.volume(), T));
    }

    print_jackknife(std::cout, results);
}
