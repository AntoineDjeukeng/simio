#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include "simio/simio.hpp"

namespace simio::analysis {

enum class Species : int { Water = 0, Na = 1, Cl = 2, SpeciesN = 3 };

inline int species_count() { return static_cast<int>(Species::SpeciesN); }

inline int species_index(Species s) { return static_cast<int>(s); }

inline int species_index_from_type(MolType t) {
    if (t == MolType::Water) return species_index(Species::Water);
    if (t == MolType::Cation) return species_index(Species::Na);
    if (t == MolType::Anion) return species_index(Species::Cl);
    return -1;
}

struct IntervalMap {
    bool inside = false;
    double u = 0.0;
    double length = 0.0;
};

inline IntervalMap map_on_pbc_interval(double xw, double mnw, double mxw, double L) {
    IntervalMap out{};
    if (L <= 0.0) return out;
    if (mnw <= mxw) {
        out.length = mxw - mnw;
        if (out.length <= 0.0) return out;
        if (xw < mnw || xw >= mxw) return out;
        out.inside = true;
        out.u = xw - mnw;
        return out;
    }
    out.length = (mxw + L) - mnw;
    if (out.length <= 0.0) return out;
    if (xw >= mnw) {
        out.inside = true;
        out.u = xw - mnw;
        return out;
    }
    if (xw < mxw) {
        out.inside = true;
        out.u = (xw + L) - mnw;
        return out;
    }
    return out;
}

inline bool in_range_pbc(double xw, double mnw, double mxw) {
    if (mnw <= mxw) return (xw >= mnw && xw < mxw);
    return (xw >= mnw) || (xw < mxw);
}

inline double interval_length_pbc(double mnw, double mxw, double L) {
    if (mnw <= mxw) return mxw - mnw;
    return (mxw + L) - mnw;
}

inline double wrapped_interval_midpoint(const Pbc3D& pbc, int axis, double mnw, double mxw, double L) {
    if (L <= 0.0) return 0.0;
    if (mnw <= mxw) return mnw + 0.5 * (mxw - mnw);
    const double u = mnw + 0.5 * ((mxw + L) - mnw);
    return pbc.wrap_pos(axis, u);
}

struct RunningStatsAll {
    std::vector<double> sum;
    std::vector<double> sumsq;

    void init(int n) {
        sum.assign(static_cast<size_t>(n), 0.0);
        sumsq.assign(static_cast<size_t>(n), 0.0);
    }

    void add(int i, double v) {
        sum[static_cast<size_t>(i)] += v;
        sumsq[static_cast<size_t>(i)] += v * v;
    }

    double mean(int i, int nframes) const {
        if (nframes <= 0) return 0.0;
        return sum[static_cast<size_t>(i)] / static_cast<double>(nframes);
    }

    double sem(int i, int nframes) const {
        if (nframes <= 1) return 0.0;
        const double s = sum[static_cast<size_t>(i)];
        const double ss = sumsq[static_cast<size_t>(i)];
        const double m = s / static_cast<double>(nframes);
        double var = (ss - static_cast<double>(nframes) * m * m) / static_cast<double>(nframes - 1);
        if (var < 0.0) var = 0.0;
        return std::sqrt(var / static_cast<double>(nframes));
    }
};

struct RunningStatsNonEmpty {
    std::vector<double> sum;
    std::vector<double> sumsq;
    std::vector<int64_t> n_used;

    void init(int n) {
        sum.assign(static_cast<size_t>(n), 0.0);
        sumsq.assign(static_cast<size_t>(n), 0.0);
        n_used.assign(static_cast<size_t>(n), 0);
    }

    void add(int i, double v) {
        const size_t idx = static_cast<size_t>(i);
        sum[idx] += v;
        sumsq[idx] += v * v;
        n_used[idx] += 1;
    }

    double mean(int i) const {
        const size_t idx = static_cast<size_t>(i);
        const int64_t n = n_used[idx];
        if (n <= 0) return 0.0;
        return sum[idx] / static_cast<double>(n);
    }

    double sem(int i) const {
        const size_t idx = static_cast<size_t>(i);
        const int64_t n = n_used[idx];
        if (n <= 1) return 0.0;
        const double s = sum[idx];
        const double ss = sumsq[idx];
        const double m = s / static_cast<double>(n);
        double var = (ss - static_cast<double>(n) * m * m) / static_cast<double>(n - 1);
        if (var < 0.0) var = 0.0;
        return std::sqrt(var / static_cast<double>(n));
    }
};

}  // namespace simio::analysis
