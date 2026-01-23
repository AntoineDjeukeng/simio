#include "mdconv/backend/writer.hpp"
#include "mdconv/diag/errors.hpp"
#include <vector>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <cmath>
#include <vector>

namespace mdconv {

void Writer::write_lammps_stub(std::ostream& os, const SystemIR& ir) const
{
    (void)ir;
    os << "LAMMPS data file via mdconv (stub)\n";
    os << "\n";
    os << "0 atoms\n";
    os << "0 atom types\n";
    os << "\n";
    os << "0.0 0.0 xlo xhi\n";
    os << "0.0 0.0 ylo yhi\n";
    os << "0.0 0.0 zlo zhi\n";
    os << "\n";
    os << "# Sections will be added in writer v1.1\n";
}

void Writer::write_lammps_stub_file(const std::string& path, const SystemIR& ir) const
{
    std::ofstream f(path.c_str(), std::ios::out | std::ios::trunc);
    if (!f) {
        throw WriteError("cannot open output file: " + path);
    }
    write_lammps_stub(f, ir);
}

void Writer::write_lammps_data(const SystemIR& ir, std::ostream& out) const
{
    const double nm_to_ang = 10.0;

    out << "LAMMPS data file via mdconv\n";
    out << "\n";
    out << ir.atoms.size() << " atoms\n";
    out << ir.bonds.size() << " bonds\n";
    out << ir.angles.size() << " angles\n";
    out << ir.dihedrals.size() << " dihedrals\n";
    out << "\n";
    out << ir.atom_types.size() << " atom types\n";
    out << ir.bond_types.size() << " bond types\n";
    out << ir.angle_types.size() << " angle types\n";
    out << ir.rb_types.size() << " dihedral types\n";
    out << "\n";

    double lx = ir.box.lx_nm;
    double ly = ir.box.ly_nm;
    double lz = ir.box.lz_nm;
    double xhi = (lx == 0.0) ? 10.0 : (lx * nm_to_ang);
    double yhi = (ly == 0.0) ? 10.0 : (ly * nm_to_ang);
    double zhi = (lz == 0.0) ? 10.0 : (lz * nm_to_ang);
    out << "0 " << xhi << " xlo xhi\n";
    out << "0 " << yhi << " ylo yhi\n";
    out << "0 " << zhi << " zlo zhi\n";
    out << "\n";

    out << "Masses\n\n";
    for (size_t i = 0; i < ir.atom_types.size(); ++i) {
        out << ir.atom_types[i].id << " " << ir.atom_types[i].mass_amu << "\n";
    }
    out << "\n";

    out << "Atoms\n\n";
    for (size_t i = 0; i < ir.atoms.size(); ++i) {
        const Atom& a = ir.atoms[i];
        out << a.id << " " << a.molecule_id << " " << a.type_id << " " << a.charge_e
            << " " << (a.x_nm * nm_to_ang)
            << " " << (a.y_nm * nm_to_ang)
            << " " << (a.z_nm * nm_to_ang) << "\n";
    }
    out << "\n";

    out << "Bonds\n\n";
    for (size_t i = 0; i < ir.bonds.size(); ++i) {
        const Bond& b = ir.bonds[i];
        out << (i + 1) << " " << b.type_id << " " << b.i << " " << b.j << "\n";
    }
    out << "\n";

    out << "Angles\n\n";
    for (size_t i = 0; i < ir.angles.size(); ++i) {
        const Angle& a = ir.angles[i];
        out << (i + 1) << " " << a.type_id << " " << a.i << " " << a.j
            << " " << a.k << "\n";
    }
    out << "\n";

    out << "Dihedrals\n\n";
    for (size_t i = 0; i < ir.dihedrals.size(); ++i) {
        const Dihedral& d = ir.dihedrals[i];
        out << (i + 1) << " " << d.type_id << " " << d.i << " " << d.j
            << " " << d.k << " " << d.l << "\n";
    }
    out << "\n";
}

void Writer::write_lammps_input_snippet(const SystemIR& ir,
                                        std::ostream& out,
                                        const Options& options) const
{
    const double kj_to_kcal = 1.0 / 4.184;
    const double nm_to_ang = 10.0;
    const double nm2_per_ang2 = 100.0;

    out << std::setprecision(17);
    out << "units real\n";
    out << "atom_style full\n";
    out << "pair_style lj/cut/coul/cut " << options.cut_ang
        << " " << options.cut_ang << "\n";
    out << "special_bonds lj 0.0 0.0 " << ir.nb_defaults.fudge_lj
        << " coul 0.0 0.0 " << ir.nb_defaults.fudge_qq << "\n";
    out << "bond_style harmonic\n";
    out << "angle_style harmonic\n";
    out << "dihedral_style ryckaert/bellemans\n";
    if (options.kspace == KspaceMode::PPPM) {
        out << "kspace_style pppm " << options.pppm_accuracy << "\n";
    }

    std::vector<double> sigma_a;
    std::vector<double> epsilon_kcal;
    sigma_a.reserve(ir.atom_types.size());
    epsilon_kcal.reserve(ir.atom_types.size());
    for (size_t i = 0; i < ir.atom_types.size(); ++i) {
        sigma_a.push_back(ir.atom_types[i].sigma_nm * nm_to_ang);
        epsilon_kcal.push_back(ir.atom_types[i].epsilon_kj * kj_to_kcal);
    }
    for (size_t i = 0; i < sigma_a.size(); ++i) {
        for (size_t j = i; j < sigma_a.size(); ++j) {
            double sigma_ij = 0.5 * (sigma_a[i] + sigma_a[j]);
            double epsilon_ij = std::sqrt(epsilon_kcal[i] * epsilon_kcal[j]);
            out << "pair_coeff " << (i + 1) << " " << (j + 1) << " "
                << epsilon_ij << " " << sigma_ij << "\n";
        }
    }

    for (size_t i = 0; i < ir.bond_types.size(); ++i) {
        const BondType& bt = ir.bond_types[i];
        double k_out = (bt.k_kj_mol_nm2 * kj_to_kcal) / nm2_per_ang2;
        double r0_out = bt.r0_nm * nm_to_ang;
        out << "bond_coeff " << bt.id << " " << k_out << " " << r0_out << "\n";
    }

    const double pi = std::acos(-1.0);
    for (size_t i = 0; i < ir.angle_types.size(); ++i) {
        const AngleType& at = ir.angle_types[i];
        double k_out = at.k_kj_mol_rad2 * kj_to_kcal;
        double theta_deg = at.theta0_rad * (180.0 / pi);
        out << "angle_coeff " << at.id << " " << k_out << " " << theta_deg << "\n";
    }

    for (size_t i = 0; i < ir.rb_types.size(); ++i) {
        const RBType& rt = ir.rb_types[i];
        out << "dihedral_coeff " << rt.id;
        for (size_t c = 0; c < rt.c_kj.size(); ++c) {
            out << " " << (rt.c_kj[c] * kj_to_kcal);
        }
        out << "\n";
    }
}

} // namespace mdconv
