// xtc_min.c
#include "xtc.h"

#include <cstdio>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

struct Setup
{
    // from topo_setup.json
    double channel_z_min = 0.902000;
    double channel_z_max = 1.802000;

    // full box lengths from gro_box / trajectory
    double box_x = 22.000000;
    double x_bin_width = 0.1;

    // filtered atom layout from topo_setup.json
    int n_cl  = 114;
    int n_na  = 138;
    int n_can = 2304;
    int n_car = 4608;
    int n_sol = 9420;

    // SOL = O H H
    int sol_atoms_per_mol = 3;

    const char* out_csv = "water_density_x.csv";

    int sol_atom_start() const
    {
        return n_cl + n_na + n_can + n_car;   // 7164
    }

    int sol_atom_end() const
    {
        return sol_atom_start() + n_sol * sol_atoms_per_mol; // 35424
    }

    double slab_thickness_z() const
    {
        return channel_z_max - channel_z_min;
    }
};

static double wrap_coord(double x, double L)
{
    double y = std::fmod(x, L);
    if (y < 0.0)
        y += L;
    return y;
}

void compute_water_density_x(const char* xtc_path)
{
    Setup set;

    const int sol_start = set.sol_atom_start();
    const int sol_end   = set.sol_atom_end();

    const int nbins = static_cast<int>(std::ceil(set.box_x / set.x_bin_width));
    if (nbins <= 0)
        throw std::runtime_error("Invalid number of x bins");

    std::vector<double> density_accum(nbins, 0.0);

    XtcTraj t = {};
    int rc = t_traj_open_auto(&t, xtc_path, 2);
    if (rc != exdrOK)
        throw std::runtime_error("Failed to open trajectory");

    int nframes = 0;

    while ((rc = xtc_read_next(&t)) == exdrOK)
    {
        const XtcFrame* fr = xtc_tail(&t);
        if (!fr)
            break;

        const double Lx = fr->box[0][0];
        const double Ly = fr->box[1][1];
        const double dz = set.slab_thickness_z();
        const double dx = set.x_bin_width;
        const double bin_volume = dx * Ly * dz;

        if (bin_volume <= 0.0)
        {
            xtc_close(&t);
            throw std::runtime_error("Invalid bin volume");
        }

        // Loop over water molecules only, take first atom of each triplet = oxygen
        for (int i = sol_start; i < sol_end; i += set.sol_atoms_per_mol)
        {
            if (i < 0 || i >= fr->natoms)
                continue;

            const double x = wrap_coord(fr->x[i][0], Lx);
            const double z = fr->x[i][2];

            if (z < set.channel_z_min || z > set.channel_z_max)
                continue;

            int bin = static_cast<int>(std::floor(x / dx));
            if (bin < 0)
                continue;
            if (bin >= nbins)
                bin = nbins - 1;

            density_accum[bin] += 1.0 / bin_volume;
        }

        ++nframes;
    }

    xtc_close(&t);

    if (rc != exdrENDOFFILE)
        throw std::runtime_error("Error while reading trajectory");

    if (nframes == 0)
        throw std::runtime_error("No frames read from trajectory");

    for (double& v : density_accum)
        v /= static_cast<double>(nframes);

    std::FILE* f = std::fopen(set.out_csv, "w");
    if (!f)
        throw std::runtime_error("Failed to open output CSV");

    std::fprintf(f, "x_center_nm,rho_water_mean\n");
    for (int b = 0; b < nbins; ++b)
    {
        const double x_lo = b * set.x_bin_width;
        const double x_hi = x_lo + set.x_bin_width;
        const double x_center = 0.5 * (x_lo + x_hi);
        std::fprintf(f, "%.6f,%.8f\n", x_center, density_accum[b]);
    }

    std::fclose(f);

    std::printf("Computed water x-density in channel z-slab\n");
    std::printf("xtc             : %s\n", xtc_path);
    std::printf("SOL atom range  : [%d, %d)\n", sol_start, sol_end);
    std::printf("water molecules : %d\n", set.n_sol);
    std::printf("frames          : %d\n", nframes);
    std::printf("z slab          : [%.3f, %.3f]\n", set.channel_z_min, set.channel_z_max);
    std::printf("dz              : %.3f\n", set.slab_thickness_z());
    std::printf("x bin width     : %.3f\n", set.x_bin_width);
    std::printf("nbins           : %d\n", nbins);
    std::printf("output          : %s\n", set.out_csv);
}

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::fprintf(stderr, "Usage: %s traj.xtc\n", argv[0]);
        return 1;
    }

    try
    {
        compute_water_density_x(argv[1]);
    }
    catch (const std::exception& e)
    {
        std::fprintf(stderr, "Error: %s\n", e.what());
        return 1;
    }

    return 0;
}
