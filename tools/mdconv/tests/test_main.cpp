#include "mdconv/gmx/parser.hpp"
#include "mdconv/validate/validator.hpp"
#include "mdconv/backend/writer.hpp"
#include "mdconv/diag/errors.hpp"
#include "mdconv/util/dump.hpp"
#include "mdconv/gmx/gro.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

static int g_failures = 0;

static std::string read_file(const std::string& path)
{
    std::ifstream in(path.c_str(), std::ios::in);
    if (!in && path.find("tests/fixtures/") == 0) {
        std::string alt = path.substr(std::string("tests/").size());
        in.clear();
        in.open(alt.c_str(), std::ios::in);
        if (!in) {
            throw std::runtime_error("failed to open file: " + path + " (also tried: " + alt + ")");
        }
    } else if (!in) {
        throw std::runtime_error("failed to open file: " + path);
    }

    std::ostringstream buffer;
    buffer << in.rdbuf();
    if (!in.good() && !in.eof()) {
        throw std::runtime_error("failed to read file: " + path);
    }
    return buffer.str();
}

static std::string replace_once(std::string s, const std::string& needle, const std::string& repl)
{
    std::string::size_type pos = s.find(needle);
    if (pos == std::string::npos) {
        throw std::runtime_error("needle not found: " + needle);
    }
    if (s.find(needle, pos + needle.size()) != std::string::npos) {
        throw std::runtime_error("needle occurs more than once: " + needle);
    }
    s.replace(pos, needle.size(), repl);
    return s;
}

static bool nearly(double a, double b, double eps)
{
    return std::fabs(a - b) <= eps;
}

static void apply_gro(mdconv::SystemIR& ir, const mdconv::GroFrame& gro)
{
    if (gro.x_nm.size() != ir.atoms.size()) {
        throw std::runtime_error("GRO atom count does not match topology");
    }
    for (size_t i = 0; i < ir.atoms.size(); ++i) {
        ir.atoms[i].x_nm = gro.x_nm[i].x_nm;
        ir.atoms[i].y_nm = gro.x_nm[i].y_nm;
        ir.atoms[i].z_nm = gro.x_nm[i].z_nm;
    }
    ir.box = gro.box_nm;
}

#define ASSERT_TRUE(cond) do { \
    if (!(cond)) { \
        std::cerr << "ASSERT_TRUE failed: " #cond \
                  << " (" << __FILE__ << ":" << __LINE__ << ")\n"; \
        ++g_failures; \
    } \
} while(0)

#define ASSERT_NO_THROW(expr) do { \
    bool thrown = false; \
    try { (void)(expr); } catch (...) { thrown = true; } \
    if (thrown) { \
        std::cerr << "ASSERT_NO_THROW failed: " #expr \
                  << " (" << __FILE__ << ":" << __LINE__ << ")\n"; \
        ++g_failures; \
    } \
} while(0)

#define ASSERT_THROWS(expr) do { \
    bool thrown = false; \
    try { (void)(expr); } catch (...) { thrown = true; } \
    if (!thrown) { \
        std::cerr << "ASSERT_THROWS failed: " #expr \
                  << " (" << __FILE__ << ":" << __LINE__ << ")\n"; \
        ++g_failures; \
    } \
} while(0)

int main()
{
    mdconv::Parser parser;
    mdconv::Validator validator;
    mdconv::Writer writer;

    // Parser throws on empty input
    ASSERT_THROWS(parser.parse_from_string(""));
    ASSERT_THROWS(parser.parse_from_string("   \n\t  "));
    ASSERT_THROWS(parser.parse_from_string("no sections here"));

    // Section scan succeeds without parsing values
    {
        mdconv::SystemIR ir = parser.parse_from_string(
            read_file("tests/fixtures/topo_valid.txt"));
        ASSERT_TRUE(ir.nb_defaults.fudge_lj == 0.5);
        ASSERT_TRUE(ir.nb_defaults.fudge_qq == 0.8333);
    }

    // Invalid combination rule throws
    {
        const std::string input = "[ defaults ]\n"
                                  "1 1 yes 0.5 0.8333\n";
        ASSERT_THROWS(parser.parse_from_string(input));
    }

    // Atomtypes parse (name, mass, sigma, epsilon)
    {
        mdconv::SystemIR ir = parser.parse_from_string(
            read_file("tests/fixtures/topo_valid.txt"));
        ASSERT_TRUE(ir.atom_types.size() == 2);
        ASSERT_TRUE(ir.atom_types[0].id == 1);
        ASSERT_TRUE(ir.atom_types[0].name == "opls_135");
        ASSERT_TRUE(ir.atom_types[0].mass_amu == 12.011);
        ASSERT_TRUE(ir.atom_types[0].sigma_nm == 0.355);
        ASSERT_TRUE(ir.atom_types[0].epsilon_kj == 0.276144);
        ASSERT_TRUE(ir.atom_types[1].id == 2);
        ASSERT_TRUE(ir.atom_types[1].name == "opls_136");
    }

    // Duplicate atomtype name throws
    {
        const std::string input = "[ defaults ]\n"
                                  "1 2 yes 0.5 0.8333\n"
                                  "[ atomtypes ]\n"
                                  "opls_135 12.011 0.0 A 0.355 0.276144\n"
                                  "opls_135 1.008 0.0 A 0.250 0.12552\n";
        ASSERT_THROWS(parser.parse_from_string(input));
    }

    // Atoms parse and validate
    {
        mdconv::SystemIR ir = parser.parse_from_string(
            read_file("tests/fixtures/topo_valid.txt"));
        mdconv::GroFrame gro = mdconv::parse_gro_from_string(
            read_file("tests/fixtures/simple.gro"));
        apply_gro(ir, gro);
        ASSERT_NO_THROW(validator.validate_or_throw(ir));
    }

    // Atoms parse with 8-field format (includes mass)
    {
        mdconv::SystemIR ir = parser.parse_from_string(
            read_file("tests/fixtures/topo_atoms_8field.txt"));
        ASSERT_TRUE(ir.atoms.size() == 4);
        ASSERT_TRUE(nearly(ir.atoms[0].charge_e, -0.120, 1e-12));
    }

    // Atom mass mismatch against [ atomtypes ] throws
    {
        std::string input = read_file("tests/fixtures/topo_atoms_8field.txt");
        input = replace_once(input,
                             "1 opls_135 1 RES C1 1 -0.120 12.011\n",
                             "1 opls_135 1 RES C1 1 -0.120 12.999\n");
        ASSERT_THROWS(parser.parse_from_string(input));
    }

    // Moleculetype + molecules assemble into a global system
    {
        std::vector<std::string> include_dirs;
        include_dirs.push_back("tests/fixtures");
        mdconv::SystemIR ir = parser.parse_from_file("tests/fixtures/system.top", include_dirs);
        ASSERT_TRUE(ir.atoms.size() == 5);
        ASSERT_TRUE(ir.bonds.size() == 2);
        ASSERT_TRUE(ir.angles.size() == 1);
        ASSERT_TRUE(ir.atoms[0].molecule_id == 1);
        ASSERT_TRUE(ir.atoms[1].molecule_id == 2);
        ASSERT_TRUE(ir.atoms[2].molecule_id == 3);
        ASSERT_TRUE(ir.atoms[4].molecule_id == 3);
    }

    // Unsupported sections fail loudly in moleculetype mode
    {
        const std::string input = "[ moleculetype ]\n"
                                  "X 1\n"
                                  "[ atoms ]\n"
                                  "1 opls_135 1 RES X 1 0.0 12.0\n"
                                  "[ nonbond_params ]\n"
                                  "X X 1 0.3 0.1\n"
                                  "[ molecules ]\n"
                                  "X 1\n";
        ASSERT_THROWS(parser.parse_from_string(input));
    }

    // Wrong atom id ordering fails validation
    {
        mdconv::SystemIR ir;
        ir.atom_types.push_back(mdconv::AtomType{1, "opls_135", 12.0, 0.3, 0.1});
        ir.atoms.push_back(mdconv::Atom{1, 1, 1, -0.1, 0.0, 0.0, 0.0});
        ir.atoms.push_back(mdconv::Atom{3, 1, 1, 0.1, 0.0, 0.0, 0.0});
        ASSERT_THROWS(validator.validate_or_throw(ir));
    }

    // Validator fails on zero atoms
    {
        mdconv::SystemIR ir;
        ASSERT_THROWS(validator.validate_or_throw(ir));
    }

    // Writer emits stub header (no throw)
    {
        mdconv::SystemIR ir;
        std::ostringstream oss;
        writer.write_lammps_stub(oss, ir);
        const std::string out = oss.str();
        ASSERT_TRUE(out.find("LAMMPS data file via mdconv (stub)") != std::string::npos);
    }

    // JSON dump contains required keys and stable header ordering
    {
        mdconv::SystemIR ir;
        ir.atom_types.push_back(mdconv::AtomType{1, "opls_135", 12.0, 0.3, 0.1});
        ir.atoms.push_back(mdconv::Atom{1, 1, 1, -0.1, 1.0, 2.0, 3.0});
        ir.box.lx_nm = 10.0;
        ir.box.ly_nm = 11.0;
        ir.box.lz_nm = 12.0;
        const std::string json = mdconv::to_json(ir);
        ASSERT_TRUE(json.rfind("{\"unit_system\":\"GMX\",\"counts\":", 0) == 0);
        ASSERT_TRUE(json.find("\"atom_types\"") != std::string::npos);
        ASSERT_TRUE(json.find("\"atoms\"") != std::string::npos);
        ASSERT_TRUE(json.find("\"nb_defaults\"") != std::string::npos);
        ASSERT_TRUE(json.find("\"box\"") != std::string::npos);
    }

    // GRO parse returns coordinates and box
    {
        mdconv::GroFrame frame = mdconv::parse_gro_from_string(
            read_file("tests/fixtures/simple.gro"));
        ASSERT_TRUE(frame.x_nm.size() == 4);
        ASSERT_TRUE(nearly(frame.x_nm[0].x_nm, 0.0, 1e-12));
        ASSERT_TRUE(nearly(frame.x_nm[1].y_nm, 0.2, 1e-12));
        ASSERT_TRUE(nearly(frame.x_nm[2].z_nm, 0.6, 1e-12));
        ASSERT_TRUE(nearly(frame.box_nm.lx_nm, 3.0, 1e-12));
        ASSERT_TRUE(nearly(frame.box_nm.ly_nm, 3.0, 1e-12));
        ASSERT_TRUE(nearly(frame.box_nm.lz_nm, 3.0, 1e-12));
    }

    // GRO coordinates and box validate when copied into IR
    {
        mdconv::SystemIR ir = parser.parse_from_string(
            read_file("tests/fixtures/topo_valid.txt"));
        mdconv::GroFrame gro = mdconv::parse_gro_from_string(
            read_file("tests/fixtures/simple.gro"));
        ASSERT_TRUE(gro.x_nm.size() == ir.atoms.size());
        for (size_t i = 0; i < ir.atoms.size(); ++i) {
            ir.atoms[i].x_nm = gro.x_nm[i].x_nm;
            ir.atoms[i].y_nm = gro.x_nm[i].y_nm;
            ir.atoms[i].z_nm = gro.x_nm[i].z_nm;
        }
        ir.box = gro.box_nm;
        ASSERT_NO_THROW(validator.validate_or_throw(ir));
    }

    // GRO zero box dimension fails validation
    {
        mdconv::SystemIR ir = parser.parse_from_string(
            read_file("tests/fixtures/topo_valid.txt"));
        std::string gro_text = read_file("tests/fixtures/simple.gro");
        gro_text = replace_once(gro_text,
                                "3.00000 3.00000 3.00000\n",
                                "0.00000 3.00000 3.00000\n");
        mdconv::GroFrame gro = mdconv::parse_gro_from_string(gro_text);
        ASSERT_TRUE(gro.x_nm.size() == ir.atoms.size());
        for (size_t i = 0; i < ir.atoms.size(); ++i) {
            ir.atoms[i].x_nm = gro.x_nm[i].x_nm;
            ir.atoms[i].y_nm = gro.x_nm[i].y_nm;
            ir.atoms[i].z_nm = gro.x_nm[i].z_nm;
        }
        ir.box = gro.box_nm;
        ASSERT_THROWS(validator.validate_or_throw(ir));
    }

    // Bondtypes and bonds parse
    {
        mdconv::SystemIR ir = parser.parse_from_string(
            read_file("tests/fixtures/topo_valid.txt"));
        mdconv::GroFrame gro = mdconv::parse_gro_from_string(
            read_file("tests/fixtures/simple.gro"));
        apply_gro(ir, gro);
        ASSERT_TRUE(ir.bond_types.size() == 1);
        ASSERT_TRUE(ir.bonds.size() == 1);
        ASSERT_NO_THROW(validator.validate_or_throw(ir));
    }

    // Unsupported bond function throws
    {
        std::string input = read_file("tests/fixtures/topo_valid.txt");
        input = replace_once(input,
                             "opls_135 opls_136 1 1000 0.1\n",
                             "opls_135 opls_136 2 1000 0.1\n");
        ASSERT_THROWS(parser.parse_from_string(input));
    }

    // Missing bondtype for pair throws
    {
        std::string input = read_file("tests/fixtures/topo_valid.txt");
        input = replace_once(input,
                             "opls_135 opls_136 1 1000 0.1\n",
                             "opls_135 opls_135 1 1000 0.1\n");
        ASSERT_THROWS(parser.parse_from_string(input));
    }

    // Duplicate bondtype pair throws
    {
        std::string input = read_file("tests/fixtures/topo_valid.txt");
        input = replace_once(input,
                             "opls_135 opls_136 1 1000 0.1\n",
                             "opls_135 opls_136 1 1000 0.1\n"
                             "opls_136 opls_135 1 900 0.2\n");
        ASSERT_THROWS(parser.parse_from_string(input));
    }

    // Duplicate bondtype pair with swapped order throws
    {
        ASSERT_THROWS(parser.parse_from_string(
            read_file("tests/fixtures/topo_dup_bondtypes_swapped.txt")));
    }

    // Invalid bond atom id throws
    {
        std::string input = read_file("tests/fixtures/topo_valid.txt");
        input = replace_once(input, "1 2 1\n", "1 5 1\n");
        ASSERT_THROWS(parser.parse_from_string(input));
    }

    // Invalid bond type id fails validation
    {
        mdconv::SystemIR ir;
        ir.atom_types.push_back(mdconv::AtomType{1, "opls_135", 12.0, 0.3, 0.1});
        ir.atoms.push_back(mdconv::Atom{1, 1, 1, -0.1, 0.0, 0.0, 0.0});
        ir.bond_types.push_back(mdconv::BondType{1, mdconv::TermForm::HarmonicBond, 1000.0, 0.1});
        ir.bonds.push_back(mdconv::Bond{1, 1, 2});
        ASSERT_THROWS(validator.validate_or_throw(ir));
    }

    // Angletypes and angles parse and validate
    {
        mdconv::SystemIR ir = parser.parse_from_string(
            read_file("tests/fixtures/topo_valid.txt"));
        mdconv::GroFrame gro = mdconv::parse_gro_from_string(
            read_file("tests/fixtures/simple.gro"));
        apply_gro(ir, gro);
        ASSERT_TRUE(ir.angle_types.size() == 1);
        ASSERT_TRUE(ir.angles.size() == 1);
        ASSERT_NO_THROW(validator.validate_or_throw(ir));
    }

    // Unsupported angletype function throws
    {
        std::string input = read_file("tests/fixtures/topo_valid.txt");
        input = replace_once(input,
                             "opls_135 opls_136 opls_135 1 100 109.5\n",
                             "opls_135 opls_136 opls_135 2 100 109.5\n");
        ASSERT_THROWS(parser.parse_from_string(input));
    }

    // Unsupported angle function throws
    {
        std::string input = read_file("tests/fixtures/topo_valid.txt");
        input = replace_once(input, "1 2 3 1\n", "1 2 3 2\n");
        ASSERT_THROWS(parser.parse_from_string(input));
    }

    // Missing angletype throws
    {
        std::string input = read_file("tests/fixtures/topo_valid.txt");
        input = replace_once(input,
                             "opls_135 opls_136 opls_135 1 100 109.5\n",
                             "opls_135 opls_135 opls_135 1 100 109.5\n");
        ASSERT_THROWS(parser.parse_from_string(input));
    }

    // Invalid angle atom id throws
    {
        std::string input = read_file("tests/fixtures/topo_valid.txt");
        input = replace_once(input, "1 2 3 1\n", "1 2 4 1\n");
        ASSERT_THROWS(parser.parse_from_string(input));
    }

    // Invalid angle type id fails validation
    {
        mdconv::SystemIR ir;
        ir.atom_types.push_back(mdconv::AtomType{1, "opls_135", 12.0, 0.3, 0.1});
        ir.atoms.push_back(mdconv::Atom{1, 1, 1, -0.1, 0.0, 0.0, 0.0});
        ir.atoms.push_back(mdconv::Atom{2, 1, 1, 0.1, 0.0, 0.0, 0.0});
        ir.atoms.push_back(mdconv::Atom{3, 1, 1, 0.1, 0.0, 0.0, 0.0});
        ir.angles.push_back(mdconv::Angle{1, 2, 3, 1});
        ASSERT_THROWS(validator.validate_or_throw(ir));
    }

    // Duplicate angletype triplet (including swapped order) throws
    {
        ASSERT_THROWS(parser.parse_from_string(
            read_file("tests/fixtures/topo_dup_angletypes.txt")));
    }

    // Dihedraltypes and dihedrals parse and validate
    {
        mdconv::SystemIR ir = parser.parse_from_string(
            read_file("tests/fixtures/topo_valid.txt"));
        mdconv::GroFrame gro = mdconv::parse_gro_from_string(
            read_file("tests/fixtures/simple.gro"));
        apply_gro(ir, gro);
        ASSERT_TRUE(ir.rb_types.size() == 1);
        ASSERT_TRUE(ir.dihedrals.size() == 1);
        ASSERT_NO_THROW(validator.validate_or_throw(ir));
    }

    // Unsupported dihedraltype function throws
    {
        std::string input = read_file("tests/fixtures/topo_valid.txt");
        input = replace_once(input,
                             "opls_135 opls_136 opls_135 opls_136 3 0 1 2 3 4 5\n",
                             "opls_135 opls_136 opls_135 opls_136 2 0 1 2 3 4 5\n");
        ASSERT_THROWS(parser.parse_from_string(input));
    }

    // Unsupported dihedral function throws
    {
        std::string input = read_file("tests/fixtures/topo_valid.txt");
        input = replace_once(input, "1 2 3 4 3\n", "1 2 3 4 1\n");
        ASSERT_THROWS(parser.parse_from_string(input));
    }

    // Missing dihedraltype throws
    {
        std::string input = read_file("tests/fixtures/topo_valid.txt");
        input = replace_once(input,
                             "opls_135 opls_136 opls_135 opls_136 3 0 1 2 3 4 5\n",
                             "opls_135 opls_135 opls_135 opls_135 3 0 1 2 3 4 5\n");
        ASSERT_THROWS(parser.parse_from_string(input));
    }

    // Invalid dihedral atom id throws
    {
        std::string input = read_file("tests/fixtures/topo_valid.txt");
        input = replace_once(input, "1 2 3 4 3\n", "1 2 3 5 3\n");
        ASSERT_THROWS(parser.parse_from_string(input));
    }

    // Duplicate dihedraltype quartet (including reverse symmetry) throws
    {
        ASSERT_THROWS(parser.parse_from_string(
            read_file("tests/fixtures/topo_dup_dihedraltypes.txt")));
    }

    // Invalid dihedral type id fails validation
    {
        mdconv::SystemIR ir;
        ir.atom_types.push_back(mdconv::AtomType{1, "opls_135", 12.0, 0.3, 0.1});
        ir.atoms.push_back(mdconv::Atom{1, 1, 1, -0.1, 0.0, 0.0, 0.0});
        ir.atoms.push_back(mdconv::Atom{2, 1, 1, 0.1, 0.0, 0.0, 0.0});
        ir.atoms.push_back(mdconv::Atom{3, 1, 1, 0.1, 0.0, 0.0, 0.0});
        ir.atoms.push_back(mdconv::Atom{4, 1, 1, 0.1, 0.0, 0.0, 0.0});
        ir.dihedrals.push_back(mdconv::Dihedral{1, 2, 3, 4, 1});
        ASSERT_THROWS(validator.validate_or_throw(ir));
    }

    // Writer data and input snippet emit expected sections and coeffs
    {
        mdconv::SystemIR ir = parser.parse_from_string(
            read_file("tests/fixtures/topo_valid.txt"));

        std::ostringstream data;
        writer.write_lammps_data(ir, data);
        const std::string data_out = data.str();
        ASSERT_TRUE(data_out.find("Masses") != std::string::npos);
        ASSERT_TRUE(data_out.find("Atoms") != std::string::npos);
        ASSERT_TRUE(data_out.find("Bonds") != std::string::npos);
        ASSERT_TRUE(data_out.find("Angles") != std::string::npos);
        ASSERT_TRUE(data_out.find("Dihedrals") != std::string::npos);
        ASSERT_TRUE(data_out.find("4 atoms") != std::string::npos);
        ASSERT_TRUE(data_out.find("1 bonds") != std::string::npos);
        ASSERT_TRUE(data_out.find("1 angles") != std::string::npos);
        ASSERT_TRUE(data_out.find("1 dihedrals") != std::string::npos);

        std::ostringstream snippet;
        mdconv::Writer::Options wopts;
        writer.write_lammps_input_snippet(ir, snippet, wopts);
        const std::string snippet_out = snippet.str();
        ASSERT_TRUE(snippet_out.find("pair_style lj/cut/coul/cut") != std::string::npos);
        ASSERT_TRUE(snippet_out.find("special_bonds lj 0.0 0.0") != std::string::npos);
        ASSERT_TRUE(snippet_out.find("pair_coeff 1 1") != std::string::npos);
        ASSERT_TRUE(snippet_out.find("pair_coeff 1 2") != std::string::npos);
        ASSERT_TRUE(snippet_out.find("pair_coeff 2 2") != std::string::npos);
        ASSERT_TRUE(snippet_out.find("bond_style harmonic") != std::string::npos);
        ASSERT_TRUE(snippet_out.find("angle_style harmonic") != std::string::npos);
        ASSERT_TRUE(snippet_out.find("dihedral_style ryckaert/bellemans") != std::string::npos);
        ASSERT_TRUE(snippet_out.find("kspace_style pppm") != std::string::npos);
        ASSERT_TRUE(snippet_out.find("bond_coeff 1") != std::string::npos);
        ASSERT_TRUE(snippet_out.find("angle_coeff 1") != std::string::npos);
        ASSERT_TRUE(snippet_out.find("dihedral_coeff 1") != std::string::npos);
    }

    // Invalid bond type id with empty bond_types fails validation
    {
        mdconv::SystemIR ir;
        ir.atom_types.push_back(mdconv::AtomType{1, "opls_135", 12.0, 0.3, 0.1});
        ir.atoms.push_back(mdconv::Atom{1, 1, 1, -0.1, 0.0, 0.0, 0.0});
        ir.atoms.push_back(mdconv::Atom{2, 1, 1, 0.1, 0.0, 0.0, 0.0});
        ir.bonds.push_back(mdconv::Bond{1, 2, 1});
        ASSERT_THROWS(validator.validate_or_throw(ir));
    }

    // Bond atom index out of range fails validation
    {
        mdconv::SystemIR ir;
        ir.atom_types.push_back(mdconv::AtomType{1, "opls_135", 12.0, 0.3, 0.1});
        ir.atoms.push_back(mdconv::Atom{1, 1, 1, -0.1, 0.0, 0.0, 0.0});
        ir.atoms.push_back(mdconv::Atom{2, 1, 1, 0.1, 0.0, 0.0, 0.0});
        ir.bond_types.push_back(mdconv::BondType{1, mdconv::TermForm::HarmonicBond, 1000.0, 0.1});
        ir.bonds.push_back(mdconv::Bond{1, 3, 1});
        ASSERT_THROWS(validator.validate_or_throw(ir));
    }

    if (g_failures == 0) {
        std::cout << "All tests passed.\n";
        return 0;
    }
    std::cerr << g_failures << " test(s) failed.\n";
    return 1;
}
