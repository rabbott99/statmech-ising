#ifndef __LATTICE_H__

#define __LATTICE_H__

#include <vector>
#include <iostream>
#include <random>

#include <cassert>

const int Nd = 2;

// Computes val + diff modulo L
inline static int cyclic_shift(int val, int diff, int L) {
    int ret = (val + diff) % L;
    // Modulo is weird in C/C++ so if val + diff < 0 then we can get a result
    // in (-L, 0) here
    if (ret < 0) {
        ret += L;
    }
    return ret;
}

template <typename T>
class Lattice {
    public:
        int L; // Size of the box
        std::vector<T> values;

        Lattice(int _L): L(_L) {
            values = std::vector<T>(volume());
        }

        int volume() const {
            int vol = 1;
            for (int mu = 0; mu < Nd; mu++) {
                vol *= L;
            }
            return vol;
        }

        std::vector<int> indexToCoord(int idx) const {
            std::vector<int> coord(Nd);
            int stride = volume();
            for (int mu = 0; mu < Nd; mu++) {
                stride /= L;
                coord[mu] = idx / stride;
                idx -= coord[mu] * stride;
            }
            return coord;
        }

        int coordToIndex(const std::vector<int> &coord) const {
            assert(coord.size() == Nd);
            int idx = 0;
            for (size_t mu = 0; mu < coord.size(); mu++) {
                if (mu > 0) {
                    idx *= L;
                }
                idx += coord[mu];
            }
            return idx;
        }

        T& get_value(const std::vector<int> &coord) {
            return values[coordToIndex(coord)];
        }

        const T& get_value_const(const std::vector<int> &coord) const {
            return values[coordToIndex(coord)];
        }

        Lattice<T> shift(int dir, int amt) const {
            assert(dir >= 0 && dir < Nd);
            Lattice<T> ret(L);
            for (int idx = 0; idx < volume(); idx++) {
                std::vector<int> coord = indexToCoord(idx);
                coord[dir] = cyclic_shift(coord[dir], amt, L);
                int out_idx = coordToIndex(coord);
                ret.values[out_idx] = values[idx];
            }
            return ret;
        }

        void pretty_print_2d(std::ostream &os) {
            for (int idx = 0; idx < volume(); idx++) {
                if (idx != 0 && idx % L == 0)
                    os << std::endl;
                os << values[idx] << " ";
            }
            os << std::endl;
        }
};


template <typename T>
Lattice<T> Cshift(const Lattice<T> &lat, int mu, int amt) {
    return lat.shift(mu, amt);
}

typedef Lattice<int> LatticeInt;
typedef Lattice<double> LatticeReal;

inline void print_spin_lat(std::ostream &os, LatticeInt &latt) {
    for (int idx = 0; idx < latt.volume(); idx++) {
        os << (latt.values[idx] > 0 ? "x" : "o");
        if (idx % latt.L == latt.L - 1) {
            os << std::endl;
        }
    }
}

LatticeInt random_lattice(int L, int seed = 0);
LatticeInt cold_lattice(int L, int val = 1);
int get_sum_neighbors(const LatticeInt &latt, int idx);

#endif
