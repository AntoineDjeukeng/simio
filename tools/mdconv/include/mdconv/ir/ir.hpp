#pragma once

#include <array>
#include <string>
#include <utility>
#include <vector>

namespace mdconv {

enum class UnitSystem { GMX };
enum class MixingRule { LorentzBerthelot };
enum class TermForm { HarmonicBond, HarmonicAngle, RBDihedral, LJ126 };

struct Vec3 {
    double x_nm{0.0};
    double y_nm{0.0};
    double z_nm{0.0};
};

struct AtomType {
    int         id = 0;          // 1..N
    std::string name;            // e.g. "opls_135"
    double      mass_amu = 0.0;
    double      sigma_nm = 0.0;
    double      epsilon_kj = 0.0;
};

struct Atom {
    int     id = 0;              // 1..Natoms
    int     molecule_id = 0;     // 1..Nmolecules
    int     type_id = 0;         // AtomType::id
    double  charge_e = 0.0;
    double  x_nm = 0.0, y_nm = 0.0, z_nm = 0.0;
};

struct BondType {
    int     id = 0;
    TermForm form = TermForm::HarmonicBond;
    double  k_kj_mol_nm2 = 0.0;
    double  r0_nm = 0.0;
};

struct Bond {
    int i = 0, j = 0;
    int type_id = 0;
};

struct AngleType {
    int     id = 0;
    TermForm form = TermForm::HarmonicAngle;
    double  k_kj_mol_rad2 = 0.0;
    double  theta0_rad = 0.0;
};

struct Angle {
    int i = 0, j = 0, k = 0;
    int type_id = 0;
};

struct RBType {
    int     id = 0;
    TermForm form = TermForm::RBDihedral;
    std::array<double, 6> c_kj{{0, 0, 0, 0, 0, 0}};
};

struct Dihedral {
    int i = 0, j = 0, k = 0, l = 0;
    int type_id = 0;
};

struct MoleculeType {
    std::string name;
    int nrexcl = 0;
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::vector<Angle> angles;
    std::vector<Dihedral> dihedrals;
    std::vector<std::pair<int, int> > exclusions;
    std::vector<std::pair<int, int> > pairs;
};

struct MoleculeCount {
    std::string name;
    int count = 0;
};

struct SystemSpec {
    std::vector<MoleculeCount> molecules;
};

struct NonbondDefaults {
    MixingRule mixing = MixingRule::LorentzBerthelot;
    double     fudge_lj = 1.0;
    double     fudge_qq = 1.0;
};

struct PairParam {
    int    type_i = 0;
    int    type_j = 0;
    double sigma_nm = 0.0;
    double epsilon_kj = 0.0;
};

struct PairScale {
    int    i = 0;
    int    j = 0;
    double scale_lj = 1.0;
    double scale_coul = 1.0;
};

struct Box {
    double lx_nm = 0.0;
    double ly_nm = 0.0;
    double lz_nm = 0.0;
};

struct SystemIR {
    UnitSystem units = UnitSystem::GMX;

    std::vector<AtomType> atom_types;
    std::vector<Atom>     atoms;

    std::vector<BondType> bond_types;
    std::vector<Bond>     bonds;

    std::vector<AngleType> angle_types;
    std::vector<Angle>     angles;

    std::vector<RBType>    rb_types;
    std::vector<Dihedral>  dihedrals;

    NonbondDefaults nb_defaults;
    std::vector<PairParam>  pair_params;
    std::vector<PairScale>  pair_scales;

    std::vector<MoleculeType> molecule_types;
    SystemSpec system;

    Box box;
};

} // namespace mdconv
