#include "mdconv/validate/validator.hpp"
#include "mdconv/diag/errors.hpp"

#include <cmath>

namespace mdconv {

void Validator::validate_or_throw(const SystemIR& ir) const
{
    if (ir.atoms.empty()) {
        throw ValidationError("system has zero atoms");
    }

    for (size_t i = 0; i < ir.atoms.size(); ++i) {
        if (ir.atoms[i].id != static_cast<int>(i + 1)) {
            throw ValidationError("atom ids must be sequential starting at 1");
        }
    }

    for (size_t i = 0; i < ir.atoms.size(); ++i) {
        int type_id = ir.atoms[i].type_id;
        if (type_id < 1 || type_id > static_cast<int>(ir.atom_types.size())) {
            throw ValidationError("atom type id out of range");
        }
        if (ir.atoms[i].molecule_id < 1) {
            throw ValidationError("atom molecule id must be positive");
        }
    }

    for (size_t i = 0; i < ir.atom_types.size(); ++i) {
        const AtomType& at = ir.atom_types[i];
        if (!std::isfinite(at.mass_amu)
            || !std::isfinite(at.sigma_nm)
            || !std::isfinite(at.epsilon_kj)) {
            throw ValidationError("atom type parameter is not finite");
        }
        if (at.mass_amu < 0.0 || at.sigma_nm < 0.0 || at.epsilon_kj < 0.0) {
            throw ValidationError("atom type parameter must be non-negative");
        }
    }

    for (size_t i = 0; i < ir.bond_types.size(); ++i) {
        const BondType& bt = ir.bond_types[i];
        if (!std::isfinite(bt.k_kj_mol_nm2) || !std::isfinite(bt.r0_nm)) {
            throw ValidationError("bond parameter is not finite");
        }
        if (bt.k_kj_mol_nm2 < 0.0 || bt.r0_nm < 0.0) {
            throw ValidationError("bond parameter must be non-negative");
        }
    }

    for (size_t i = 0; i < ir.angle_types.size(); ++i) {
        const AngleType& at = ir.angle_types[i];
        if (!std::isfinite(at.k_kj_mol_rad2) || !std::isfinite(at.theta0_rad)) {
            throw ValidationError("angle parameter is not finite");
        }
        if (at.k_kj_mol_rad2 < 0.0 || at.theta0_rad < 0.0) {
            throw ValidationError("angle parameter must be non-negative");
        }
    }

    for (size_t i = 0; i < ir.rb_types.size(); ++i) {
        const RBType& rt = ir.rb_types[i];
        for (size_t c = 0; c < rt.c_kj.size(); ++c) {
            if (!std::isfinite(rt.c_kj[c])) {
                throw ValidationError("dihedral parameter is not finite");
            }
        }
    }

    for (size_t i = 0; i < ir.atoms.size(); ++i) {
        const Atom& atom = ir.atoms[i];
        if (!std::isfinite(atom.x_nm)
            || !std::isfinite(atom.y_nm)
            || !std::isfinite(atom.z_nm)) {
            throw ValidationError("atom coordinate is not finite");
        }
    }

    if (!std::isfinite(ir.box.lx_nm)
        || !std::isfinite(ir.box.ly_nm)
        || !std::isfinite(ir.box.lz_nm)) {
        throw ValidationError("box is not finite");
    }
    if (ir.box.lx_nm <= 0.0 || ir.box.ly_nm <= 0.0 || ir.box.lz_nm <= 0.0) {
        throw ValidationError("box dimensions must be positive");
    }

    for (size_t i = 0; i < ir.bonds.size(); ++i) {
        int atom_i = ir.bonds[i].i;
        int atom_j = ir.bonds[i].j;
        if (atom_i < 1 || atom_j < 1
            || atom_i > static_cast<int>(ir.atoms.size())
            || atom_j > static_cast<int>(ir.atoms.size())) {
            throw ValidationError("bond atom index out of range");
        }

        int type_id = ir.bonds[i].type_id;
        if (type_id < 1 || type_id > static_cast<int>(ir.bond_types.size())) {
            throw ValidationError("bond type id out of range");
        }
    }

    for (size_t i = 0; i < ir.angles.size(); ++i) {
        int atom_i = ir.angles[i].i;
        int atom_j = ir.angles[i].j;
        int atom_k = ir.angles[i].k;
        if (atom_i < 1 || atom_j < 1 || atom_k < 1
            || atom_i > static_cast<int>(ir.atoms.size())
            || atom_j > static_cast<int>(ir.atoms.size())
            || atom_k > static_cast<int>(ir.atoms.size())) {
            throw ValidationError("angle atom index out of range");
        }

        int type_id = ir.angles[i].type_id;
        if (type_id < 1 || type_id > static_cast<int>(ir.angle_types.size())) {
            throw ValidationError("angle type id out of range");
        }
    }

    for (size_t i = 0; i < ir.dihedrals.size(); ++i) {
        int atom_i = ir.dihedrals[i].i;
        int atom_j = ir.dihedrals[i].j;
        int atom_k = ir.dihedrals[i].k;
        int atom_l = ir.dihedrals[i].l;
        if (atom_i < 1 || atom_j < 1 || atom_k < 1 || atom_l < 1
            || atom_i > static_cast<int>(ir.atoms.size())
            || atom_j > static_cast<int>(ir.atoms.size())
            || atom_k > static_cast<int>(ir.atoms.size())
            || atom_l > static_cast<int>(ir.atoms.size())) {
            throw ValidationError("dihedral atom index out of range");
        }

        int type_id = ir.dihedrals[i].type_id;
        if (type_id < 1 || type_id > static_cast<int>(ir.rb_types.size())) {
            throw ValidationError("dihedral type id out of range");
        }
    }
}

} // namespace mdconv
