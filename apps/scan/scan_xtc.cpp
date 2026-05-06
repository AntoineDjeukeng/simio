#include "xtc.h"
#include "json_utils.hpp"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include <string>
#include <cstdint>

struct Setup
{
    int ion_index_0 = 0;   // inclusive, 0-based
    int ion_index_1 = 0;   // exclusive, 0-based

    double channel_x_min = 0.0;
    double channel_x_max = 0.0;

    double channel_z_min = 0.0;
    double channel_z_max = 0.0;

    double x_pad = 1.0;
    double x_bin_width = 0.1;
    double tmax_ps = 20000.0;

    int n_groups = 4;

    Setup(int i0, int i1,
          double cx_min, double cx_max,
          double cz_min, double cz_max)
        : ion_index_0(i0),
          ion_index_1(i1),
          channel_x_min(cx_min),
          channel_x_max(cx_max),
          channel_z_min(cz_min),
          channel_z_max(cz_max)
    {}

    Setup() = default;

    double roi_x_min() const { return channel_x_min - x_pad; }
    double roi_x_max() const { return channel_x_max + x_pad; }
    double z_target() const { return 0.5 * (channel_z_min + channel_z_max); }
};

Setup make_setup(const json_utils::json& j, const std::string& mol_name)
{
    json_utils::IntRange ion = json_utils::find_atom_range(j, mol_name);

    double cx_min = json_utils::find_value<double>(j, {"channel", "min", "0"});
    double cx_max = json_utils::find_value<double>(j, {"channel", "max", "0"});

    double cz_min = json_utils::find_value<double>(j, {"channel", "min", "2"});
    double cz_max = json_utils::find_value<double>(j, {"channel", "max", "2"});

    return Setup(
        ion.min,
        ion.max + 1,   // convert inclusive max to exclusive end
        cx_min,
        cx_max,
        cz_min,
        cz_max
    );
}

struct BestPerBin
{
    bool found = false;

    int bin_index = -1;
    double x_lo = 0.0;
    double x_hi = 0.0;

    int frame_id = -1;
    int step = -1;
    double time = 0.0;

    int atom_index = -1;   // 0-based
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    double dz_abs = std::numeric_limits<double>::infinity();
};

static void update_best_per_bin(
    std::vector<BestPerBin>& bins,
    int bin_index,
    const Setup& set,
    int frame_id,
    int step,
    double time,
    int atom_index,
    double x,
    double y,
    double z,
    double dz_abs)
{
    BestPerBin& best = bins[bin_index];
    const double eps = 1e-9;

    const bool better =
        !best.found ||
        (dz_abs < best.dz_abs - eps) ||
        (std::fabs(dz_abs - best.dz_abs) <= eps && frame_id > best.frame_id);

    if (better)
    {
        best.found = true;
        best.bin_index = bin_index;
        best.x_lo = set.roi_x_min() + bin_index * set.x_bin_width;
        best.x_hi = best.x_lo + set.x_bin_width;

        best.frame_id = frame_id;
        best.step = step;
        best.time = time;

        best.atom_index = atom_index;
        best.x = x;
        best.y = y;
        best.z = z;
        best.dz_abs = dz_abs;
    }
}

static void print_selected_candidates(
    const std::vector<BestPerBin>& best_bins,
    int nbins,
    int n_groups)
{
    std::printf("\nSelected %d candidates (range-covered, latest per segment):\n", n_groups);

    for (int g = 0; g < n_groups; ++g)
    {
        const int b0 = (g * nbins) / n_groups;
        const int b1 = ((g + 1) * nbins) / n_groups;

        int best_idx = -1;

        for (int b = b0; b < b1; ++b)
        {
            if (!best_bins[b].found)
                continue;

            if (best_idx < 0 ||
                best_bins[b].frame_id > best_bins[best_idx].frame_id ||
                (best_bins[b].frame_id == best_bins[best_idx].frame_id &&
                 best_bins[b].dz_abs < best_bins[best_idx].dz_abs))
            {
                best_idx = b;
            }
        }

        if (best_idx < 0)
        {
            std::printf("group %d: no candidate found\n", g + 1);
            continue;
        }

        const BestPerBin& c = best_bins[best_idx];

        std::printf(
            "group %d  bin [%5.3f, %5.3f]  "
            "frame=%d step=%d time=%g  "
            "atom=%d(0-based) %d(1-based)  "
            "x=%7.3f y=%7.3f z=%7.3f  |z-z0|=%8.5f\n",
            g + 1,
            c.x_lo, c.x_hi,
            c.frame_id, c.step, c.time,
            c.atom_index, c.atom_index + 1,
            c.x, c.y, c.z, c.dz_abs
        );
    }
}

void print_xtc_content(const char* path, const Setup& set, const char* mol_name)
{
    XtcTraj t = {};
    int rc = t_traj_open_auto(&t, path, 2);
    if (rc != exdrOK) {
        throw std::runtime_error("Failed to open trajectory");
    }

    const double roi_min = set.roi_x_min();
    const double roi_max = set.roi_x_max();
    const double z0 = set.z_target();

    if (roi_max <= roi_min) {
        xtc_close(&t);
        throw std::runtime_error("Invalid ROI in x");
    }

    if (set.ion_index_1 <= set.ion_index_0) {
        xtc_close(&t);
        throw std::runtime_error("Invalid ion index range");
    }

    if (set.n_groups <= 0) {
        xtc_close(&t);
        throw std::runtime_error("Invalid number of groups");
    }

    const double roi_width = roi_max - roi_min;
    const int nbins = static_cast<int>(std::ceil(roi_width / set.x_bin_width));

    if (nbins <= 0) {
        xtc_close(&t);
        throw std::runtime_error("Invalid number of bins");
    }

    std::vector<BestPerBin> best_bins(nbins);

    int frame_id = 0;
    bool stopped_at_target_time = false;

    while ((rc = xtc_read_next(&t)) == exdrOK)
    {
        const XtcFrame* fr = xtc_tail(&t);
        if (!fr)
            break;

        if (fr->natoms < set.ion_index_1) {
            xtc_close(&t);
            throw std::runtime_error("Frame does not contain enough atoms for requested ion block");
        }

        for (int i = set.ion_index_0; i < set.ion_index_1; ++i)
        {
            const double x = fr->x[i][0];
            const double y = fr->x[i][1];
            const double z = fr->x[i][2];

            if (x < roi_min || x > roi_max)
                continue;

            const double dz_abs = std::fabs(z - z0);

            int bin_index = static_cast<int>(std::floor((x - roi_min) / set.x_bin_width));

            if (bin_index < 0)
                continue;
            if (bin_index >= nbins)
                bin_index = nbins - 1;

            update_best_per_bin(
                best_bins,
                bin_index,
                set,
                frame_id,
                fr->step,
                fr->time,
                i,
                x,
                y,
                z,
                dz_abs
            );
        }

        if (fr->time >= set.tmax_ps)
        {
            std::printf("Reached %.3f ps at frame %d (step %d)\n",
                        set.tmax_ps, frame_id, fr->step);
            stopped_at_target_time = true;
            break;
        }

        ++frame_id;
    }

    std::printf("Ion scan near target z plane\n");
    std::printf("file          : %s\n", path);
    std::printf("species       : %s\n", mol_name);
    std::printf("ion block     : [%d, %d)  (0-based)\n", set.ion_index_0, set.ion_index_1);
    std::printf("ion block     : [%d, %d]  (1-based)\n", set.ion_index_0 + 1, set.ion_index_1);
    std::printf("channel_x     : [%g, %g]\n", set.channel_x_min, set.channel_x_max);
    std::printf("channel_z     : [%g, %g]\n", set.channel_z_min, set.channel_z_max);
    std::printf("z_target      : %.6f\n", z0);
    std::printf("roi_x         : [%g, %g]\n", roi_min, roi_max);
    std::printf("x_bin_width   : %g\n", set.x_bin_width);
    std::printf("nbins         : %d\n", nbins);
    std::printf("groups        : %d\n", set.n_groups);

    int nfound = 0;
    for (const auto& b : best_bins)
    {
        if (!b.found)
            continue;

        ++nfound;
        std::printf(
            "bin [%5.3f, %5.3f]  "
            "frame=%d step=%d time=%g  "
            "atom=%d(0-based) %d(1-based)  "
            "x=%7.3f y=%7.3f z=%7.3f  |z-z0|=%8.5f\n",
            b.x_lo, b.x_hi,
            b.frame_id, b.step, b.time,
            b.atom_index, b.atom_index + 1,
            b.x, b.y, b.z, b.dz_abs
        );
    }

    std::printf("filled bins   : %d / %d\n", nfound, nbins);

    print_selected_candidates(best_bins, nbins, set.n_groups);

    xtc_close(&t);

    if (!stopped_at_target_time && rc != exdrENDOFFILE) {
        throw std::runtime_error("Error while reading trajectory");
    }
}

int main(int argc, char** argv)
{
    if (argc < 4)
    {
        std::fprintf(stderr,
            "Usage: %s file.xtc topo_setup.json MOLE_NAME\n",
            argv[0]);
        return 1;
    }

    try {
        const char* xtc_file = argv[1];
        const char* json_file = argv[2];
        const char* mol_name = argv[3];

        json_utils::json j = json_utils::load_json(json_file);
        Setup set = make_setup(j, mol_name);

        print_xtc_content(xtc_file, set, mol_name);
    }
    catch (const std::exception& e) {
        std::fprintf(stderr, "Error: %s\n", e.what());
        return 1;
    }

    return 0;
}
