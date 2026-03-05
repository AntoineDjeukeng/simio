// xtc_min.c
#include "xtc.h"

#include "xtc.h"
#include <cstdio>
#include <stdexcept>
#include <vector>

#include <string>
#include <vector>
#include <cstdint>

using AtomIndex = int32_t;

struct Molecule {
    std::vector<AtomIndex> atoms;
};

struct MoleculeType 
{
    std::string name;
    std::vector<Molecule> molecules;
};

struct Topology 
{
    std::vector<MoleculeType> types;
};




void print_xtc_content(const char* path) 
{
    XtcTraj t = {};
    int rc = t_traj_open_auto(&t, path, 2); // small window is enough
    if (rc != exdrOK) {
        throw std::runtime_error("Failed to open trajectory");
    }

    int frame_id = 0;
    while ((rc = xtc_read_next(&t)) == exdrOK) {
        const XtcFrame* fr = xtc_tail(&t);
        if (!fr) break;

        std::printf("frame=%d step=%d time=%g natoms=%d\n",
                    frame_id++, fr->step, fr->time, fr->natoms);

        for (int i = 0; i < fr->natoms; ++i) {
            std::printf("  atom %d: %g %g %g\n",
                        i, fr->x[i][0], fr->x[i][1], fr->x[i][2]);
            if(i>5)
                break;
        }
    }

    xtc_close(&t);

    if (rc != exdrENDOFFILE) {
        throw std::runtime_error("Error while reading trajectory");
    }
}

int main(int argc, char** argv) {
    if (argc < 2) return 1;
    print_xtc_content(argv[1]);
    return 0;
}
