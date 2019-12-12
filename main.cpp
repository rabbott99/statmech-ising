#include <iostream>

#include "Lattice.h"
#include "evolution.h"
#include "measure.h"

template <typename T>
T average(const std::vector<T> &values) {
    T total = values[0];
    for (size_t i = 1; i < values.size(); i++)
        total += values[i];
    return total / values.size();
}

template <typename T>
std::vector<T> jackknife_resample(const std::vector<T> &samples) {
    assert(samples.size() > 1);
    T total = samples[0];
    for (size_t i = 1; i < samples.size(); i++) total += samples[i];

    std::vector<T> ret(samples.size(), total);
    for (size_t i = 0; i < samples.size(); i++) {
        ret[i] -= samples[i];
        ret[i] *= 1.0 / (samples.size() - 1);
    }

    return ret;
}

void print_result(std::ostream &os, const PhysicalResults &res) {
    os << res.energy << " " << res.heat_capacity << " "
        << res.susceptibility;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: ./ising L" << std::endl;
        exit(1);
    }
    const int L = atoi(argv[1]);

    const int thermalization_iterations = 1e5;
    const int measurement_separation = 10;
    const int num_measurements = 3e4;

    const LatticeInt start = random_lattice(L);


    const double Tmin = 0.015;
    const double Tmax = 4.5;
    const double Tstep = 0.015;
    int num_temps = std::lround((Tmax - Tmin) / Tstep) + 1;

    std::vector<PhysicalResults> output, errors;
    output.resize(num_temps);
    errors.resize(num_temps);
#pragma omp parallel for
    for (int i = 0; i < num_temps; i++) {
        double T = Tmin + i * Tstep;
        double beta = 1.0 / T;
        LatticeInt latt = start;
        std::vector<MeasurementResults> measurements = measurement_run(latt,
                beta, 0.0, thermalization_iterations,
                measurement_separation, num_measurements);

        std::vector<MeasurementResults> jackknife = jackknife_resample(measurements);
        std::vector<PhysicalResults> results;
        results.reserve(jackknife.size());
        for (auto sample: jackknife) {
            results.push_back(physical_results(sample, start.volume(), T));
        }

        output[i] = average(results);

        for (auto res: results) {
            PhysicalResults tmp = res - output[i];
            tmp = tmp * tmp;
            errors[i] += tmp;
        }
        errors[i] *= (errors.size() - 1) / errors.size();
        errors[i] = errors[i].sqrt();
    }

    for (int i = 0; i < num_temps; i++) {
        double T = Tmin + i * Tstep;
        std::cout << T << " ";
        print_result(std::cout, output[i]);
        std::cout << " ";
        print_result(std::cout, errors[i]);
        std::cout << std::endl;
    }
}
