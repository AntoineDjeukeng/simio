#pragma once

#include "simio/simio.hpp"

namespace simio::intrinsics {

struct WaterDipoleResult {
    Vec3d ref_wrapped{};
    Vec3d com_wrapped{};
    Vec3d dipole{};
    bool has_dipole = false;
};

inline WaterDipoleResult compute_water_dipole(const Frame& fr, int atom_o, int atom_h1, int atom_h2) {
    WaterDipoleResult out;

    Vec3d rO = fr.atoms.pos(atom_o);
    Vec3d rH1 = fr.atoms.pos(atom_h1);
    Vec3d rH2 = fr.atoms.pos(atom_h2);
    fr.pbc.wrap_pos3(rO);

    const Vec3d h1u = pbc_unwrap_to_ref(fr.pbc, rO, rH1);
    const Vec3d h2u = pbc_unwrap_to_ref(fr.pbc, rO, rH2);
    const Vec3d com_u = water_com_pbc_unwrapped(fr.pbc, rO, rH1, rH2);

    Vec3d com_w = com_u;
    fr.pbc.wrap_pos3(com_w);

    out.ref_wrapped = rO;
    out.com_wrapped = com_w;

    if (!fr.atoms.q.empty()) {
        auto q = [&](int ai) -> double { return fr.atoms.q[(size_t)ai]; };
        Vec3d mu{};
        mu = mu + (rO - com_u) * q(atom_o);
        mu = mu + (h1u - com_u) * q(atom_h1);
        mu = mu + (h2u - com_u) * q(atom_h2);
        out.dipole = mu;
        out.has_dipole = true;
    }

    return out;
}

}  // namespace simio::intrinsics

