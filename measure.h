#ifndef __MEASURE_H__

#define __MEASURE_H__

#include "Lattice.h"

class MeasurementResults {
    public:
        double energy, energy_squared;
        double magnetization, magnetization_squared;

        void operator +=(const MeasurementResults &other) {
            energy += other.energy;
            energy_squared += other.energy_squared;
            magnetization += other.magnetization;
            magnetization_squared += other.magnetization_squared;
        }

        void operator -=(const MeasurementResults &other) {
            energy -= other.energy;
            energy_squared -= other.energy_squared;
            magnetization -= other.magnetization;
            magnetization_squared -= other.magnetization_squared;
        }

        MeasurementResults operator +(const MeasurementResults &other) const {
            MeasurementResults ret = *this;
            ret += other;
            return ret;
        }

        void operator *=(const double &c) {
            energy *= c;
            energy_squared *= c;
            magnetization *= c;
            magnetization_squared *= c;
        }

        MeasurementResults operator *(const double &c) const {
            MeasurementResults ret = *this;
            ret *= c;
            return ret;
        }

        void operator /=(const double &c) {
            energy /= c;
            energy_squared /= c;
            magnetization /= c;
            magnetization_squared /= c;
        }

        MeasurementResults operator /(const double &c) const {
            MeasurementResults ret = *this;
            ret /= c;
            return ret;
        }

        MeasurementResults operator -() const {
            MeasurementResults ret = *this;
            ret *= -1.0;
            return ret;
        }

        MeasurementResults operator -(const MeasurementResults &other) const {
            MeasurementResults ret = *this;
            ret -= other;
            return ret;
        }
};

inline std::ostream& operator <<(std::ostream &os, const MeasurementResults& m) {
    os << "E = " << m.energy << std::endl;
    os << "E^2 = " << m.energy_squared << std::endl;
    os << "M = " << m.magnetization << std::endl;
    os << "M^2 = " << m.magnetization_squared << std::endl;
    return os;
}

struct PhysicalResults {
    double energy;
    double heat_capacity;
    double magnetization;
    double susceptibility;

    void operator +=(const PhysicalResults &other) {
        energy += other.energy;
        heat_capacity += other.heat_capacity;
        susceptibility += other.susceptibility;
        magnetization += other.magnetization;
    }

    PhysicalResults operator +(const PhysicalResults &other) const {
        PhysicalResults ret = *this;
        ret += other;
        return ret;
    }

    void operator -=(const PhysicalResults &other) {
        energy -= other.energy;
        heat_capacity -= other.heat_capacity;
        susceptibility -= other.susceptibility;
        magnetization -= other.magnetization;
    }

    PhysicalResults operator -(const PhysicalResults &other) const {
        PhysicalResults ret = *this;
        ret -= other;
        return ret;
    }

    void operator /=(const double &c) {
        energy /= c;
        heat_capacity /= c;
        susceptibility /= c;
        magnetization /= c;
    }

    PhysicalResults operator /(const double &c) const {
        PhysicalResults ret = *this;
        ret /= c;
        return ret;
    }

    void operator *=(const double &c) {
        energy *= c;
        heat_capacity *= c;
        susceptibility *= c;
        magnetization *= c;
    }

    PhysicalResults operator *(const double &c) const {
        PhysicalResults ret = *this;
        ret *= c;
        return ret;
    }

    void operator *=(const PhysicalResults &other) {
        energy *= other.energy;
        heat_capacity *= other.heat_capacity;
        susceptibility *= other.susceptibility;
        magnetization *= other.magnetization;
    }

    PhysicalResults operator *(const PhysicalResults &other) const {
        PhysicalResults ret = *this;
        ret *= other;
        return ret;
    }

    PhysicalResults sqrt() const {
        PhysicalResults ret = *this;
        ret.energy = std::sqrt(ret.energy);
        ret.heat_capacity = std::sqrt(ret.heat_capacity);
        ret.susceptibility = std::sqrt(ret.susceptibility);
        ret.magnetization = std::sqrt(ret.magnetization);
        return ret;
    }
};

inline std::ostream& operator <<(std::ostream &os, const PhysicalResults& r) {
    os << "E = " << r.energy << std::endl;
    os << "C = " << r.heat_capacity << std::endl;
    os << "magnetization = " << r.magnetization << std::endl;
    os << "chi = " << r.susceptibility << std::endl;
    return os;
}

PhysicalResults physical_results(const MeasurementResults &measurements,
        int lattice_volume, double T);

MeasurementResults measure(const LatticeInt &latt);

std::vector<MeasurementResults> measurement_run(LatticeInt &latt,
        double beta, double muB,
        int thermalization_trajectories, int measurement_separation,
        int num_measurements);

#endif
