#include "mdconv/gmx/parser.hpp"
#include "mdconv/diag/errors.hpp"
#include "mdconv/gmx/expander.hpp"

#include <array>
#include <cctype>
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace mdconv {

struct SectionLine {
    std::string section;
    std::string line;
    size_t index;
};

bool Parser::is_all_whitespace(const std::string& s)
{
    for (size_t i = 0; i < s.size(); ++i) {
        if (!std::isspace(static_cast<unsigned char>(s[i])))
            return false;
    }
    return true;
}

static bool parse_int(const std::string& s, int& out)
{
    std::istringstream in(s);
    in >> out;
    return in && in.eof();
}

static bool parse_double(const std::string& s, double& out)
{
    std::istringstream in(s);
    in >> out;
    return in && in.eof();
}

static std::string strip_gmx_comment(const std::string& line)
{
    std::string::size_type pos = line.find(';');
    if (pos == std::string::npos) {
        return line;
    }
    return line.substr(0, pos);
}

static bool has_tokens(const std::string& line)
{
    std::istringstream in(line);
    std::string tok;
    return static_cast<bool>(in >> tok);
}

static std::string trim_whitespace(const std::string& line)
{
    size_t start = 0;
    while (start < line.size() && std::isspace(static_cast<unsigned char>(line[start]))) {
        ++start;
    }
    size_t end = line.size();
    while (end > start && std::isspace(static_cast<unsigned char>(line[end - 1]))) {
        --end;
    }
    return line.substr(start, end - start);
}

static bool parse_section_header(const std::string& trimmed, std::string& out_section)
{
    if (trimmed.size() < 3 || trimmed.front() != '[' || trimmed.back() != ']') {
        return false;
    }
    size_t name_start = 1;
    while (name_start + 1 < trimmed.size()
           && std::isspace(static_cast<unsigned char>(trimmed[name_start]))) {
        ++name_start;
    }
    size_t name_end = trimmed.size() - 1;
    while (name_end > name_start
           && std::isspace(static_cast<unsigned char>(trimmed[name_end - 1]))) {
        --name_end;
    }
    out_section = trimmed.substr(name_start, name_end - name_start);
    return !out_section.empty();
}

static std::string expand_topology(const std::string& root,
                                   const std::vector<std::string>& include_dirs,
                                   bool dump_to_stdout)
{
    mdconv::Expander ex(root);
    for (size_t i = 0; i < include_dirs.size(); ++i) {
        ex.add_include_dir(include_dirs[i]);
    }

    std::ostringstream out;
    std::string line;
    while (ex.getline(line)) {
        if (dump_to_stdout) {
            std::cout << line << "\n";
        }
        out << line << "\n";
    }
    return out.str();
}

static std::string format_section_error(const std::string& section,
                                        size_t index,
                                        const std::string& line,
                                        const std::string& msg)
{
    std::ostringstream out;
    out << msg << " in [" << section << "] at line " << index << ": " << line;
    return out.str();
}

static void ensure_supported_section(const std::string& section)
{
    if (section == "defaults"
        || section == "atomtypes"
        || section == "bondtypes"
        || section == "angletypes"
        || section == "dihedraltypes"
        || section == "moleculetype"
        || section == "atoms"
        || section == "bonds"
        || section == "angles"
        || section == "dihedrals"
        || section == "exclusions"
        || section == "pairs"
        || section == "molecules"
        || section == "system") {
        return;
    }
    throw ParseError("unsupported section: [ " + section + " ]");
}

static void parse_molecule_blocks(const std::string& input,
                                  const std::map<std::string, int>& atomtype_by_name,
                                  const std::vector<AtomType>& atom_types,
                                  const std::map<std::pair<int, int>, int>& bondtype_by_pair,
                                  const std::map<std::tuple<int, int, int>, int>& angletype_by_triplet,
                                  const std::map<std::tuple<int, int, int, int>, int>& dihedraltype_by_quad,
                                  std::vector<MoleculeType>& molecule_types,
                                  SystemSpec& system_spec)
{
    std::istringstream in(input);
    std::string line;
    std::string current_section;
    MoleculeType* current_mol = nullptr;
    bool expect_moleculetype_header = false;

    while (std::getline(in, line)) {
        std::string raw = strip_gmx_comment(line);
        std::string trimmed = trim_whitespace(raw);
        if (trimmed.empty()) {
            continue;
        }
        if (!trimmed.empty() && trimmed[0] == '#') {
            throw ParseError("preprocessor directives are not supported: " + trimmed);
        }

        std::string header;
        if (parse_section_header(trimmed, header)) {
            current_section = header;
            ensure_supported_section(current_section);
            if (current_section == "moleculetype") {
                expect_moleculetype_header = true;
                current_mol = nullptr;
            } else if (current_section == "molecules" || current_section == "system") {
                current_mol = nullptr;
            }
            continue;
        }

        if (current_section == "defaults"
            || current_section == "atomtypes"
            || current_section == "bondtypes"
            || current_section == "angletypes"
            || current_section == "dihedraltypes") {
            continue;
        }

        if (current_section == "moleculetype") {
            if (!expect_moleculetype_header) {
                throw ParseError("unexpected line in [ moleculetype ]: " + trimmed);
            }
            std::istringstream header_line(trimmed);
            std::vector<std::string> tokens;
            std::string tok;
            while (header_line >> tok) {
                tokens.push_back(tok);
            }
            if (tokens.size() < 2) {
                throw ParseError("invalid [ moleculetype ] header: " + trimmed);
            }
            MoleculeType mt;
            mt.name = tokens[0];
            if (!parse_int(tokens[1], mt.nrexcl)) {
                throw ParseError("invalid [ moleculetype ] header: " + trimmed);
            }
            molecule_types.push_back(mt);
            current_mol = &molecule_types.back();
            expect_moleculetype_header = false;
            continue;
        }

        if (current_section == "molecules") {
            std::istringstream mol_line(trimmed);
            std::vector<std::string> tokens;
            std::string tok;
            while (mol_line >> tok) {
                tokens.push_back(tok);
            }
            if (tokens.size() != 2) {
                throw ParseError("invalid [ molecules ] entry: " + trimmed);
            }
            int count = 0;
            if (!parse_int(tokens[1], count) || count < 0) {
                throw ParseError("invalid [ molecules ] entry: " + trimmed);
            }
            MoleculeCount mc;
            mc.name = tokens[0];
            mc.count = count;
            system_spec.molecules.push_back(mc);
            continue;
        }

        if (current_section == "system") {
            continue;
        }

        if (current_mol == nullptr) {
            throw ParseError("section [ " + current_section + " ] appears outside [ moleculetype ]");
        }

        if (current_section == "atoms") {
            std::istringstream atom_line(trimmed);
            std::vector<std::string> tokens;
            std::string field;
            while (atom_line >> field) {
                tokens.push_back(field);
            }
            if (tokens.size() != 7 && tokens.size() != 8) {
                throw ParseError("unsupported [ atoms ] format (expected 7 or 8 fields): " + trimmed);
            }

            int atom_id = 0;
            int resnr = 0;
            int cgnr = 0;
            double charge = 0.0;
            if (!parse_int(tokens[0], atom_id)
                || !parse_int(tokens[2], resnr)
                || !parse_int(tokens[5], cgnr)
                || !parse_double(tokens[6], charge)) {
                throw ParseError("invalid [ atoms ] entry: " + trimmed);
            }
            if (atom_id != static_cast<int>(current_mol->atoms.size() + 1)) {
                throw ParseError("atom ids must be sequential starting at 1 in [ atoms ]: " + trimmed);
            }

            const std::string& type_name = tokens[1];
            std::map<std::string, int>::const_iterator type_it =
                atomtype_by_name.find(type_name);
            if (type_it == atomtype_by_name.end()) {
                throw ParseError("unknown atomtype in [ atoms ]: " + type_name);
            }
            int type_id = type_it->second;

            if (tokens.size() == 8) {
                double mass = 0.0;
                if (!parse_double(tokens[7], mass)) {
                    throw ParseError("invalid [ atoms ] entry: " + trimmed);
                }
                if (type_id < 1 || type_id > static_cast<int>(atom_types.size())) {
                    throw ParseError("atom type id out of range in [ atoms ]");
                }
                const AtomType& at = atom_types[type_id - 1];
                if (std::fabs(at.mass_amu - mass) > 1e-4) {
                    throw ParseError("atom mass disagrees with [ atomtypes ] for type " + type_name);
                }
            }

            Atom atom;
            atom.id = atom_id;
            atom.type_id = type_id;
            atom.charge_e = charge;
            atom.molecule_id = 0;
            current_mol->atoms.push_back(atom);
            continue;
        }

        if (current_section == "bonds") {
            std::istringstream bond_line(trimmed);
            std::vector<std::string> tokens;
            std::string field;
            while (bond_line >> field) {
                tokens.push_back(field);
            }
            if (tokens.size() < 3) {
                throw ParseError("invalid [ bonds ] entry: " + trimmed);
            }
            if (tokens[2] != "1") {
                throw ParseError("unsupported bond function in [ bonds ]: " + trimmed);
            }
            int ai = 0;
            int aj = 0;
            if (!parse_int(tokens[0], ai) || !parse_int(tokens[1], aj)) {
                throw ParseError("invalid [ bonds ] entry: " + trimmed);
            }
            if (ai < 1 || aj < 1
                || ai > static_cast<int>(current_mol->atoms.size())
                || aj > static_cast<int>(current_mol->atoms.size())) {
                throw ParseError("invalid [ bonds ] atom index: " + trimmed);
            }

            int type_i = current_mol->atoms[ai - 1].type_id;
            int type_j = current_mol->atoms[aj - 1].type_id;
            int first = type_i < type_j ? type_i : type_j;
            int second = type_i < type_j ? type_j : type_i;
            std::map<std::pair<int, int>, int>::const_iterator bondtype_it =
                bondtype_by_pair.find(std::make_pair(first, second));
            if (bondtype_it == bondtype_by_pair.end()) {
                throw ParseError("no matching [ bondtypes ] for [ bonds ]: " + trimmed);
            }

            Bond bond;
            bond.i = ai;
            bond.j = aj;
            bond.type_id = bondtype_it->second;
            current_mol->bonds.push_back(bond);
            continue;
        }

        if (current_section == "angles") {
            std::istringstream angle_line(trimmed);
            std::vector<std::string> tokens;
            std::string field;
            while (angle_line >> field) {
                tokens.push_back(field);
            }
            if (tokens.size() < 4) {
                throw ParseError("invalid [ angles ] entry: " + trimmed);
            }
            if (tokens[3] != "1") {
                throw ParseError("unsupported angle function in [ angles ]: " + trimmed);
            }
            int ai = 0;
            int aj = 0;
            int ak = 0;
            if (!parse_int(tokens[0], ai)
                || !parse_int(tokens[1], aj)
                || !parse_int(tokens[2], ak)) {
                throw ParseError("invalid [ angles ] entry: " + trimmed);
            }
            if (ai < 1 || aj < 1 || ak < 1
                || ai > static_cast<int>(current_mol->atoms.size())
                || aj > static_cast<int>(current_mol->atoms.size())
                || ak > static_cast<int>(current_mol->atoms.size())) {
                throw ParseError("invalid [ angles ] atom index: " + trimmed);
            }

            int type_i = current_mol->atoms[ai - 1].type_id;
            int type_j = current_mol->atoms[aj - 1].type_id;
            int type_k = current_mol->atoms[ak - 1].type_id;
            std::map<std::tuple<int, int, int>, int>::const_iterator angle_it =
                angletype_by_triplet.find(std::make_tuple(type_i, type_j, type_k));
            if (angle_it == angletype_by_triplet.end()) {
                throw ParseError("no matching [ angletypes ] for [ angles ]: " + trimmed);
            }

            Angle angle;
            angle.i = ai;
            angle.j = aj;
            angle.k = ak;
            angle.type_id = angle_it->second;
            current_mol->angles.push_back(angle);
            continue;
        }

        if (current_section == "dihedrals") {
            std::istringstream dihedral_line(trimmed);
            std::vector<std::string> tokens;
            std::string field;
            while (dihedral_line >> field) {
                tokens.push_back(field);
            }
            if (tokens.size() < 5) {
                throw ParseError("invalid [ dihedrals ] entry: " + trimmed);
            }
            if (tokens[4] != "3") {
                throw ParseError("unsupported dihedral function in [ dihedrals ]: " + trimmed);
            }
            int ai = 0;
            int aj = 0;
            int ak = 0;
            int al = 0;
            if (!parse_int(tokens[0], ai)
                || !parse_int(tokens[1], aj)
                || !parse_int(tokens[2], ak)
                || !parse_int(tokens[3], al)) {
                throw ParseError("invalid [ dihedrals ] entry: " + trimmed);
            }
            if (ai < 1 || aj < 1 || ak < 1 || al < 1
                || ai > static_cast<int>(current_mol->atoms.size())
                || aj > static_cast<int>(current_mol->atoms.size())
                || ak > static_cast<int>(current_mol->atoms.size())
                || al > static_cast<int>(current_mol->atoms.size())) {
                throw ParseError("invalid [ dihedrals ] atom index: " + trimmed);
            }

            int type_i = current_mol->atoms[ai - 1].type_id;
            int type_j = current_mol->atoms[aj - 1].type_id;
            int type_k = current_mol->atoms[ak - 1].type_id;
            int type_l = current_mol->atoms[al - 1].type_id;
            std::map<std::tuple<int, int, int, int>, int>::const_iterator dihedral_it =
                dihedraltype_by_quad.find(std::make_tuple(type_i, type_j, type_k, type_l));
            if (dihedral_it == dihedraltype_by_quad.end()) {
                throw ParseError("no matching [ dihedraltypes ] for [ dihedrals ]: " + trimmed);
            }

            Dihedral dihedral;
            dihedral.i = ai;
            dihedral.j = aj;
            dihedral.k = ak;
            dihedral.l = al;
            dihedral.type_id = dihedral_it->second;
            current_mol->dihedrals.push_back(dihedral);
            continue;
        }

        if (current_section == "exclusions") {
            std::istringstream excl_line(trimmed);
            std::vector<std::string> tokens;
            std::string field;
            while (excl_line >> field) {
                tokens.push_back(field);
            }
            if (tokens.size() < 2) {
                throw ParseError("invalid [ exclusions ] entry: " + trimmed);
            }
            int ai = 0;
            if (!parse_int(tokens[0], ai)) {
                throw ParseError("invalid [ exclusions ] entry: " + trimmed);
            }
            for (size_t t = 1; t < tokens.size(); ++t) {
                int aj = 0;
                if (!parse_int(tokens[t], aj)) {
                    throw ParseError("invalid [ exclusions ] entry: " + trimmed);
                }
                current_mol->exclusions.push_back(std::make_pair(ai, aj));
            }
            continue;
        }

        if (current_section == "pairs") {
            std::istringstream pair_line(trimmed);
            std::vector<std::string> tokens;
            std::string field;
            while (pair_line >> field) {
                tokens.push_back(field);
            }
            if (tokens.size() < 2) {
                throw ParseError("invalid [ pairs ] entry: " + trimmed);
            }
            int ai = 0;
            int aj = 0;
            if (!parse_int(tokens[0], ai) || !parse_int(tokens[1], aj)) {
                throw ParseError("invalid [ pairs ] entry: " + trimmed);
            }
            current_mol->pairs.push_back(std::make_pair(ai, aj));
            continue;
        }
    }
}

static SystemIR assemble_system(const SystemIR& template_ir,
                                const std::vector<MoleculeType>& molecule_types,
                                const SystemSpec& system_spec)
{
    SystemIR out = template_ir;
    out.atoms.clear();
    out.bonds.clear();
    out.angles.clear();
    out.dihedrals.clear();
    out.molecule_types = molecule_types;
    out.system = system_spec;

    std::map<std::string, const MoleculeType*> mol_by_name;
    for (size_t i = 0; i < molecule_types.size(); ++i) {
        mol_by_name[molecule_types[i].name] = &molecule_types[i];
    }

    int next_atom_id = 1;
    int next_molecule_id = 1;
    for (size_t i = 0; i < system_spec.molecules.size(); ++i) {
        const MoleculeCount& mc = system_spec.molecules[i];
        std::map<std::string, const MoleculeType*>::const_iterator it =
            mol_by_name.find(mc.name);
        if (it == mol_by_name.end()) {
            throw ParseError("unknown moleculetype in [ molecules ]: " + mc.name);
        }
        const MoleculeType* mt = it->second;
        for (int count = 0; count < mc.count; ++count) {
            int atom_offset = next_atom_id - 1;
            for (size_t a = 0; a < mt->atoms.size(); ++a) {
                Atom atom = mt->atoms[a];
                atom.id = next_atom_id++;
                atom.molecule_id = next_molecule_id;
                out.atoms.push_back(atom);
            }
            for (size_t b = 0; b < mt->bonds.size(); ++b) {
                Bond bond = mt->bonds[b];
                bond.i += atom_offset;
                bond.j += atom_offset;
                out.bonds.push_back(bond);
            }
            for (size_t a = 0; a < mt->angles.size(); ++a) {
                Angle angle = mt->angles[a];
                angle.i += atom_offset;
                angle.j += atom_offset;
                angle.k += atom_offset;
                out.angles.push_back(angle);
            }
            for (size_t d = 0; d < mt->dihedrals.size(); ++d) {
                Dihedral dihedral = mt->dihedrals[d];
                dihedral.i += atom_offset;
                dihedral.j += atom_offset;
                dihedral.k += atom_offset;
                dihedral.l += atom_offset;
                out.dihedrals.push_back(dihedral);
            }
            ++next_molecule_id;
        }
    }
    return out;
}

SystemIR Parser::parse_from_string(const std::string& input) const
{
    if (input.empty() || is_all_whitespace(input)) {
        throw ParseError("empty input");
    }

    std::map<std::string, std::vector<SectionLine> > sections;
    std::istringstream in(input);
    std::string line;
    std::string current_section;
    bool has_moleculetype = false;
    bool has_molecules = false;

    while (std::getline(in, line)) {
        size_t comment_pos = line.find(';');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        size_t start = 0;
        while (start < line.size() && std::isspace(static_cast<unsigned char>(line[start]))) {
            ++start;
        }
        size_t end = line.size();
        while (end > start && std::isspace(static_cast<unsigned char>(line[end - 1]))) {
            --end;
        }
        std::string trimmed = line.substr(start, end - start);
        if (trimmed.empty()) {
            continue;
        }

        std::string header;
        if (parse_section_header(trimmed, header)) {
            current_section = header;
            if (!current_section.empty()) {
                sections[current_section];
            }
            continue;
        }

        if (!current_section.empty()) {
            size_t idx = sections[current_section].size();
            SectionLine entry;
            entry.section = current_section;
            entry.line = trimmed;
            entry.index = idx;
            sections[current_section].push_back(entry);
        }
    }

    if (sections.empty()) {
        throw ParseError("no sections found");
    }
    has_moleculetype = sections.find("moleculetype") != sections.end();
    has_molecules = sections.find("molecules") != sections.end();

    // Stub v1: build an empty IR but with canonical defaults.
    // Real parsing will populate atom_types, atoms, bonds, etc.
    SystemIR ir;
    std::map<std::string, std::vector<SectionLine> >::const_iterator it =
        sections.find("defaults");
    if (it == sections.end()) {
        throw ParseError("missing [ defaults ]");
    }
    if (it->second.empty()) {
        throw ParseError("invalid [ defaults ]");
    }

    const SectionLine& defaults_entry = it->second.front();
    std::string defaults_raw = strip_gmx_comment(defaults_entry.line);
    std::istringstream defaults_line(defaults_raw);
    std::vector<std::string> tokens;
    std::string tok;
    while (defaults_line >> tok) {
        tokens.push_back(tok);
    }
    if (tokens.size() < 4) {
        throw ParseError(format_section_error("defaults",
                                              defaults_entry.index,
                                              defaults_entry.line,
                                              "invalid [ defaults ]"));
    }

    int comb_rule = 0;
    if (!parse_int(tokens[1], comb_rule)) {
        throw ParseError(format_section_error("defaults",
                                              defaults_entry.index,
                                              defaults_entry.line,
                                              "invalid [ defaults ]"));
    }
    if (comb_rule != 2) {
        throw ParseError(format_section_error("defaults",
                                              defaults_entry.index,
                                              defaults_entry.line,
                                              "unsupported combination rule"));
    }

    double fudge_lj = 0.0;
    double fudge_qq = 0.0;
    if (!parse_double(tokens[tokens.size() - 2], fudge_lj)
        || !parse_double(tokens[tokens.size() - 1], fudge_qq)) {
        throw ParseError(format_section_error("defaults",
                                              defaults_entry.index,
                                              defaults_entry.line,
                                              "invalid [ defaults ]"));
    }

    ir.nb_defaults.mixing = MixingRule::LorentzBerthelot;
    ir.nb_defaults.fudge_lj = fudge_lj;
    ir.nb_defaults.fudge_qq = fudge_qq;

    it = sections.find("atomtypes");
    if (it == sections.end()) {
        throw ParseError("missing [ atomtypes ]");
    }
    int next_atomtype_id = 1;
    std::map<std::string, bool> seen_names;
    std::map<std::string, int> atomtype_by_name;
    for (size_t i = 0; i < it->second.size(); ++i) {
        const SectionLine& atomtype_entry = it->second[i];
        std::string raw = strip_gmx_comment(atomtype_entry.line);
        if (!has_tokens(raw)) {
            continue;
        }
        std::istringstream atomtype_line(raw);
        std::vector<std::string> tokens;
        std::string field;
        while (atomtype_line >> field) {
            tokens.push_back(field);
        }
        if (tokens.size() != 6 && tokens.size() != 8) {
            throw ParseError(format_section_error("atomtypes",
                                                  atomtype_entry.index,
                                                  atomtype_entry.line,
                                                  "unsupported [ atomtypes ] format (expected 6 or 8 fields)"));
        }

        std::string name;
        std::string ptype;
        double mass = 0.0;
        double sigma = 0.0;
        double epsilon = 0.0;

        if (tokens.size() == 6) {
            name = tokens[0];
            ptype = tokens[3];
            if (!parse_double(tokens[1], mass)
                || !parse_double(tokens[4], sigma)
                || !parse_double(tokens[5], epsilon)) {
                throw ParseError(format_section_error("atomtypes",
                                                      atomtype_entry.index,
                                                      atomtype_entry.line,
                                                      "invalid numeric field in [ atomtypes ]"));
            }
        } else {
            name = tokens[0];
            ptype = tokens[5];
            if (!parse_double(tokens[3], mass)
                || !parse_double(tokens[6], sigma)
                || !parse_double(tokens[7], epsilon)) {
                throw ParseError(format_section_error("atomtypes",
                                                      atomtype_entry.index,
                                                      atomtype_entry.line,
                                                      "invalid numeric field in [ atomtypes ]"));
            }
        }

        if (ptype != "A") {
            throw ParseError(format_section_error("atomtypes",
                                                  atomtype_entry.index,
                                                  atomtype_entry.line,
                                                  "unsupported ptype in [ atomtypes ] (only A/12-6 supported)"));
        }

        if (seen_names.find(name) != seen_names.end()) {
            throw ParseError(format_section_error("atomtypes",
                                                  atomtype_entry.index,
                                                  atomtype_entry.line,
                                                  "duplicate atomtype: " + name));
        }
        seen_names[name] = true;

        if (sigma < 0.0 || epsilon < 0.0) {
            throw ParseError(format_section_error("atomtypes",
                                                  atomtype_entry.index,
                                                  atomtype_entry.line,
                                                  "invalid [ atomtypes ]"));
        }

        AtomType at;
        at.id = next_atomtype_id++;
        at.name = name;
        at.mass_amu = mass;
        at.sigma_nm = sigma;
        at.epsilon_kj = epsilon;
        ir.atom_types.push_back(at);
        atomtype_by_name[name] = at.id;
    }

    std::map<std::pair<int, int>, int> bondtype_by_pair;
    std::map<std::tuple<int, int, int>, int> angletype_by_triplet;
    std::map<std::tuple<int, int, int, int>, int> dihedraltype_by_quad;

    it = sections.find("bondtypes");
    if (it == sections.end()) {
        throw ParseError("missing [ bondtypes ]");
    }
    int next_bondtype_id = 1;
    for (size_t i = 0; i < it->second.size(); ++i) {
        const SectionLine& bondtype_entry = it->second[i];
        std::string raw = strip_gmx_comment(bondtype_entry.line);
        if (!has_tokens(raw)) {
            continue;
        }
        std::istringstream bondtype_line(raw);
        std::vector<std::string> tokens;
        std::string field;
        while (bondtype_line >> field) {
            tokens.push_back(field);
        }
        if (tokens.size() < 5) {
            throw ParseError(format_section_error("bondtypes",
                                                  bondtype_entry.index,
                                                  bondtype_entry.line,
                                                  "invalid [ bondtypes ]"));
        }

        std::map<std::string, int>::const_iterator type_i_it =
            atomtype_by_name.find(tokens[0]);
        std::map<std::string, int>::const_iterator type_j_it =
            atomtype_by_name.find(tokens[1]);
        if (type_i_it == atomtype_by_name.end() || type_j_it == atomtype_by_name.end()) {
            throw ParseError(format_section_error("bondtypes",
                                                  bondtype_entry.index,
                                                  bondtype_entry.line,
                                                  "invalid [ bondtypes ]"));
        }

        if (tokens[2] != "1") {
            throw ParseError(format_section_error("bondtypes",
                                                  bondtype_entry.index,
                                                  bondtype_entry.line,
                                                  "unsupported bond function"));
        }

        double k = 0.0;
        double r0 = 0.0;
        if (!parse_double(tokens[3], k) || !parse_double(tokens[4], r0)) {
            throw ParseError(format_section_error("bondtypes",
                                                  bondtype_entry.index,
                                                  bondtype_entry.line,
                                                  "invalid [ bondtypes ]"));
        }

        BondType bt;
        bt.id = next_bondtype_id++;
        bt.form = TermForm::HarmonicBond;
        bt.k_kj_mol_nm2 = k;
        bt.r0_nm = r0;
        ir.bond_types.push_back(bt);

        int type_i = type_i_it->second;
        int type_j = type_j_it->second;
        int first = type_i < type_j ? type_i : type_j;
        int second = type_i < type_j ? type_j : type_i;
        std::pair<int, int> key = std::make_pair(first, second);
        if (bondtype_by_pair.find(key) != bondtype_by_pair.end()) {
            throw ParseError(format_section_error("bondtypes",
                                                  bondtype_entry.index,
                                                  bondtype_entry.line,
                                                  "duplicate [ bondtypes ] for pair: "
                                                      + tokens[0] + "-" + tokens[1]));
        }
        bondtype_by_pair[key] = bt.id;
    }

    it = sections.find("angletypes");
    if (it == sections.end()) {
        throw ParseError("missing [ angletypes ]");
    }
    int next_angletype_id = 1;
    for (size_t i = 0; i < it->second.size(); ++i) {
        const SectionLine& angletype_entry = it->second[i];
        std::string raw = strip_gmx_comment(angletype_entry.line);
        if (!has_tokens(raw)) {
            continue;
        }
        std::istringstream angletype_line(raw);
        std::vector<std::string> tokens;
        std::string field;
        while (angletype_line >> field) {
            tokens.push_back(field);
        }
        if (tokens.size() < 6) {
            throw ParseError(format_section_error("angletypes",
                                                  angletype_entry.index,
                                                  angletype_entry.line,
                                                  "invalid [ angletypes ]"));
        }

        std::map<std::string, int>::const_iterator type_i_it =
            atomtype_by_name.find(tokens[0]);
        std::map<std::string, int>::const_iterator type_j_it =
            atomtype_by_name.find(tokens[1]);
        std::map<std::string, int>::const_iterator type_k_it =
            atomtype_by_name.find(tokens[2]);
        if (type_i_it == atomtype_by_name.end()
            || type_j_it == atomtype_by_name.end()
            || type_k_it == atomtype_by_name.end()) {
            throw ParseError(format_section_error("angletypes",
                                                  angletype_entry.index,
                                                  angletype_entry.line,
                                                  "invalid [ angletypes ]"));
        }

        if (tokens[3] != "1") {
            throw ParseError(format_section_error("angletypes",
                                                  angletype_entry.index,
                                                  angletype_entry.line,
                                                  "unsupported angle function"));
        }

        double k = 0.0;
        double theta_deg = 0.0;
        if (!parse_double(tokens[4], k) || !parse_double(tokens[5], theta_deg)) {
            throw ParseError(format_section_error("angletypes",
                                                  angletype_entry.index,
                                                  angletype_entry.line,
                                                  "invalid [ angletypes ]"));
        }

        const double pi = std::acos(-1.0);
        // v1 assumes input angles are degrees and converts to radians.
        double theta_rad = theta_deg * (pi / 180.0);

        AngleType at;
        at.id = next_angletype_id++;
        at.form = TermForm::HarmonicAngle;
        at.k_kj_mol_rad2 = k;
        at.theta0_rad = theta_rad;
        ir.angle_types.push_back(at);

        int type_i = type_i_it->second;
        int type_j = type_j_it->second;
        int type_k = type_k_it->second;
        std::tuple<int, int, int> forward_key = std::make_tuple(type_i, type_j, type_k);
        std::tuple<int, int, int> reverse_key = std::make_tuple(type_k, type_j, type_i);
        if (angletype_by_triplet.find(forward_key) != angletype_by_triplet.end()
            || angletype_by_triplet.find(reverse_key) != angletype_by_triplet.end()) {
            throw ParseError(format_section_error("angletypes",
                                                  angletype_entry.index,
                                                  angletype_entry.line,
                                                  "duplicate [ angletypes ] for triplet: "
                                                      + tokens[0] + "-" + tokens[1] + "-" + tokens[2]));
        }
        angletype_by_triplet[forward_key] = at.id;
        angletype_by_triplet[reverse_key] = at.id;
    }

    it = sections.find("dihedraltypes");
    if (it == sections.end()) {
        throw ParseError("missing [ dihedraltypes ]");
    }
    int next_rbtype_id = 1;
    for (size_t i = 0; i < it->second.size(); ++i) {
        const SectionLine& dihedraltype_entry = it->second[i];
        std::string raw = strip_gmx_comment(dihedraltype_entry.line);
        if (!has_tokens(raw)) {
            continue;
        }
        std::istringstream dihedraltype_line(raw);
        std::vector<std::string> tokens;
        std::string field;
        while (dihedraltype_line >> field) {
            tokens.push_back(field);
        }
        if (tokens.size() < 11) {
            throw ParseError(format_section_error("dihedraltypes",
                                                  dihedraltype_entry.index,
                                                  dihedraltype_entry.line,
                                                  "invalid [ dihedraltypes ]"));
        }

        std::map<std::string, int>::const_iterator type_i_it =
            atomtype_by_name.find(tokens[0]);
        std::map<std::string, int>::const_iterator type_j_it =
            atomtype_by_name.find(tokens[1]);
        std::map<std::string, int>::const_iterator type_k_it =
            atomtype_by_name.find(tokens[2]);
        std::map<std::string, int>::const_iterator type_l_it =
            atomtype_by_name.find(tokens[3]);
        if (type_i_it == atomtype_by_name.end()
            || type_j_it == atomtype_by_name.end()
            || type_k_it == atomtype_by_name.end()
            || type_l_it == atomtype_by_name.end()) {
            throw ParseError(format_section_error("dihedraltypes",
                                                  dihedraltype_entry.index,
                                                  dihedraltype_entry.line,
                                                  "invalid [ dihedraltypes ]"));
        }

        if (tokens[4] != "3") {
            throw ParseError(format_section_error("dihedraltypes",
                                                  dihedraltype_entry.index,
                                                  dihedraltype_entry.line,
                                                  "unsupported dihedral function"));
        }

        std::array<double, 6> c_kj;
        for (size_t c = 0; c < 6; ++c) {
            if (!parse_double(tokens[5 + c], c_kj[c])) {
                throw ParseError(format_section_error("dihedraltypes",
                                                      dihedraltype_entry.index,
                                                      dihedraltype_entry.line,
                                                      "invalid [ dihedraltypes ]"));
            }
        }

        RBType rb;
        rb.id = next_rbtype_id++;
        rb.form = TermForm::RBDihedral;
        rb.c_kj = c_kj;
        ir.rb_types.push_back(rb);

        int type_i = type_i_it->second;
        int type_j = type_j_it->second;
        int type_k = type_k_it->second;
        int type_l = type_l_it->second;
        std::tuple<int, int, int, int> forward_key =
            std::make_tuple(type_i, type_j, type_k, type_l);
        std::tuple<int, int, int, int> reverse_key =
            std::make_tuple(type_l, type_k, type_j, type_i);
        if (dihedraltype_by_quad.find(forward_key) != dihedraltype_by_quad.end()
            || dihedraltype_by_quad.find(reverse_key) != dihedraltype_by_quad.end()) {
            throw ParseError(format_section_error("dihedraltypes",
                                                  dihedraltype_entry.index,
                                                  dihedraltype_entry.line,
                                                  "duplicate [ dihedraltypes ] for quartet: "
                                                      + tokens[0] + "-" + tokens[1] + "-"
                                                      + tokens[2] + "-" + tokens[3]));
        }
        dihedraltype_by_quad[forward_key] = rb.id;
        dihedraltype_by_quad[reverse_key] = rb.id;
    }

    if (!has_moleculetype && !has_molecules) {
        it = sections.find("atoms");
        if (it == sections.end()) {
            throw ParseError("missing [ atoms ]");
        }
        int next_atom_id = 1;
        for (size_t i = 0; i < it->second.size(); ++i) {
            const SectionLine& atom_entry = it->second[i];
            std::string raw = strip_gmx_comment(atom_entry.line);
            if (!has_tokens(raw)) {
                continue;
            }
            std::istringstream atom_line(raw);
            std::vector<std::string> tokens;
            std::string field;
            while (atom_line >> field) {
                tokens.push_back(field);
            }
            if (tokens.size() != 7 && tokens.size() != 8) {
                throw ParseError(format_section_error("atoms",
                                                      atom_entry.index,
                                                      atom_entry.line,
                                                      "unsupported [ atoms ] format (expected 7 or 8 fields)"));
            }

            int atom_id = 0;
            int resnr = 0;
            int cgnr = 0;
            double charge = 0.0;
            if (!parse_int(tokens[0], atom_id)
                || !parse_int(tokens[2], resnr)
                || !parse_int(tokens[5], cgnr)
                || !parse_double(tokens[6], charge)) {
                throw ParseError(format_section_error("atoms",
                                                      atom_entry.index,
                                                      atom_entry.line,
                                                      "invalid [ atoms ]"));
            }
            if (atom_id != next_atom_id) {
                throw ParseError(format_section_error("atoms",
                                                      atom_entry.index,
                                                      atom_entry.line,
                                                      "invalid [ atoms ]"));
            }

            const std::string& type_name = tokens[1];
            std::map<std::string, int>::const_iterator type_it =
                atomtype_by_name.find(type_name);
            if (type_it == atomtype_by_name.end()) {
                throw ParseError(format_section_error("atoms",
                                                      atom_entry.index,
                                                      atom_entry.line,
                                                      "unknown atomtype in [ atoms ]: " + type_name));
            }
            int type_id = type_it->second;

            if (tokens.size() == 8) {
                double mass = 0.0;
                if (!parse_double(tokens[7], mass)) {
                    throw ParseError(format_section_error("atoms",
                                                          atom_entry.index,
                                                          atom_entry.line,
                                                          "invalid [ atoms ]"));
                }
                const AtomType& at = ir.atom_types[type_id - 1];
                if (std::fabs(at.mass_amu - mass) > 1e-4) {
                    throw ParseError(format_section_error("atoms",
                                                          atom_entry.index,
                                                          atom_entry.line,
                                                          "atom mass disagrees with [ atomtypes ] for type "
                                                              + type_name));
                }
            }

            Atom atom;
            atom.id = atom_id;
            atom.type_id = type_id;
            atom.charge_e = charge;
            atom.molecule_id = 1;
            atom.x_nm = 0.0;
            atom.y_nm = 0.0;
            atom.z_nm = 0.0;
            ir.atoms.push_back(atom);
            ++next_atom_id;
        }

        it = sections.find("bonds");
        if (it == sections.end()) {
            throw ParseError("missing [ bonds ]");
        }
        for (size_t i = 0; i < it->second.size(); ++i) {
            const SectionLine& bond_entry = it->second[i];
            std::string raw = strip_gmx_comment(bond_entry.line);
            if (!has_tokens(raw)) {
                continue;
            }
            std::istringstream bond_line(raw);
            std::vector<std::string> tokens;
            std::string field;
            while (bond_line >> field) {
                tokens.push_back(field);
            }
            if (tokens.size() < 3) {
                throw ParseError(format_section_error("bonds",
                                                      bond_entry.index,
                                                      bond_entry.line,
                                                      "invalid [ bonds ]"));
            }
            if (tokens[2] != "1") {
                throw ParseError(format_section_error("bonds",
                                                      bond_entry.index,
                                                      bond_entry.line,
                                                      "unsupported bond function"));
            }

            int ai = 0;
            int aj = 0;
            if (!parse_int(tokens[0], ai) || !parse_int(tokens[1], aj)) {
                throw ParseError(format_section_error("bonds",
                                                      bond_entry.index,
                                                      bond_entry.line,
                                                      "invalid [ bonds ]"));
            }
            if (ai < 1 || aj < 1
                || ai > static_cast<int>(ir.atoms.size())
                || aj > static_cast<int>(ir.atoms.size())) {
                throw ParseError(format_section_error("bonds",
                                                      bond_entry.index,
                                                      bond_entry.line,
                                                      "invalid [ bonds ]"));
            }

            int type_i = ir.atoms[ai - 1].type_id;
            int type_j = ir.atoms[aj - 1].type_id;
            int first = type_i < type_j ? type_i : type_j;
            int second = type_i < type_j ? type_j : type_i;
            std::map<std::pair<int, int>, int>::const_iterator bondtype_it =
                bondtype_by_pair.find(std::make_pair(first, second));
            if (bondtype_it == bondtype_by_pair.end()) {
                throw ParseError(format_section_error("bonds",
                                                      bond_entry.index,
                                                      bond_entry.line,
                                                      "no matching [ bondtypes ]"));
            }

            Bond bond;
            bond.i = ai;
            bond.j = aj;
            bond.type_id = bondtype_it->second;
            ir.bonds.push_back(bond);
        }

        it = sections.find("angles");
        if (it == sections.end()) {
            throw ParseError("missing [ angles ]");
        }
        for (size_t i = 0; i < it->second.size(); ++i) {
            const SectionLine& angle_entry = it->second[i];
            std::string raw = strip_gmx_comment(angle_entry.line);
            if (!has_tokens(raw)) {
                continue;
            }
            std::istringstream angle_line(raw);
            std::vector<std::string> tokens;
            std::string field;
            while (angle_line >> field) {
                tokens.push_back(field);
            }
            if (tokens.size() < 4) {
                throw ParseError(format_section_error("angles",
                                                      angle_entry.index,
                                                      angle_entry.line,
                                                      "invalid [ angles ]"));
            }
            if (tokens[3] != "1") {
                throw ParseError(format_section_error("angles",
                                                      angle_entry.index,
                                                      angle_entry.line,
                                                      "unsupported angle function"));
            }

            int ai = 0;
            int aj = 0;
            int ak = 0;
            if (!parse_int(tokens[0], ai)
                || !parse_int(tokens[1], aj)
                || !parse_int(tokens[2], ak)) {
                throw ParseError(format_section_error("angles",
                                                      angle_entry.index,
                                                      angle_entry.line,
                                                      "invalid [ angles ]"));
            }
            if (ai < 1 || aj < 1 || ak < 1
                || ai > static_cast<int>(ir.atoms.size())
                || aj > static_cast<int>(ir.atoms.size())
                || ak > static_cast<int>(ir.atoms.size())) {
                throw ParseError(format_section_error("angles",
                                                      angle_entry.index,
                                                      angle_entry.line,
                                                      "invalid [ angles ]"));
            }

            int type_i = ir.atoms[ai - 1].type_id;
            int type_j = ir.atoms[aj - 1].type_id;
            int type_k = ir.atoms[ak - 1].type_id;
            std::map<std::tuple<int, int, int>, int>::const_iterator angle_it =
                angletype_by_triplet.find(std::make_tuple(type_i, type_j, type_k));
            if (angle_it == angletype_by_triplet.end()) {
                throw ParseError(format_section_error("angles",
                                                      angle_entry.index,
                                                      angle_entry.line,
                                                      "no matching [ angletypes ]"));
            }

            Angle angle;
            angle.i = ai;
            angle.j = aj;
            angle.k = ak;
            angle.type_id = angle_it->second;
            ir.angles.push_back(angle);
        }

        it = sections.find("dihedrals");
        if (it == sections.end()) {
            throw ParseError("missing [ dihedrals ]");
        }
        for (size_t i = 0; i < it->second.size(); ++i) {
            const SectionLine& dihedral_entry = it->second[i];
            std::string raw = strip_gmx_comment(dihedral_entry.line);
            if (!has_tokens(raw)) {
                continue;
            }
            std::istringstream dihedral_line(raw);
            std::vector<std::string> tokens;
            std::string field;
            while (dihedral_line >> field) {
                tokens.push_back(field);
            }
            if (tokens.size() < 5) {
                throw ParseError(format_section_error("dihedrals",
                                                      dihedral_entry.index,
                                                      dihedral_entry.line,
                                                      "invalid [ dihedrals ]"));
            }
            if (tokens[4] != "3") {
                throw ParseError(format_section_error("dihedrals",
                                                      dihedral_entry.index,
                                                      dihedral_entry.line,
                                                      "unsupported dihedral function"));
            }

            int ai = 0;
            int aj = 0;
            int ak = 0;
            int al = 0;
            if (!parse_int(tokens[0], ai)
                || !parse_int(tokens[1], aj)
                || !parse_int(tokens[2], ak)
                || !parse_int(tokens[3], al)) {
                throw ParseError(format_section_error("dihedrals",
                                                      dihedral_entry.index,
                                                      dihedral_entry.line,
                                                      "invalid [ dihedrals ]"));
            }
            if (ai < 1 || aj < 1 || ak < 1 || al < 1
                || ai > static_cast<int>(ir.atoms.size())
                || aj > static_cast<int>(ir.atoms.size())
                || ak > static_cast<int>(ir.atoms.size())
                || al > static_cast<int>(ir.atoms.size())) {
                throw ParseError(format_section_error("dihedrals",
                                                      dihedral_entry.index,
                                                      dihedral_entry.line,
                                                      "invalid [ dihedrals ]"));
            }

            int type_i = ir.atoms[ai - 1].type_id;
            int type_j = ir.atoms[aj - 1].type_id;
            int type_k = ir.atoms[ak - 1].type_id;
            int type_l = ir.atoms[al - 1].type_id;
            std::map<std::tuple<int, int, int, int>, int>::const_iterator dihedral_it =
                dihedraltype_by_quad.find(std::make_tuple(type_i, type_j, type_k, type_l));
            if (dihedral_it == dihedraltype_by_quad.end()) {
                throw ParseError(format_section_error("dihedrals",
                                                      dihedral_entry.index,
                                                      dihedral_entry.line,
                                                      "no matching [ dihedraltypes ]"));
            }

            Dihedral dihedral;
            dihedral.i = ai;
            dihedral.j = aj;
            dihedral.k = ak;
            dihedral.l = al;
            dihedral.type_id = dihedral_it->second;
            ir.dihedrals.push_back(dihedral);
        }
        return ir;
    }

    std::vector<MoleculeType> molecule_types;
    SystemSpec system_spec;
    parse_molecule_blocks(input,
                          atomtype_by_name,
                          ir.atom_types,
                          bondtype_by_pair,
                          angletype_by_triplet,
                          dihedraltype_by_quad,
                          molecule_types,
                          system_spec);
    if (system_spec.molecules.empty()) {
        throw ParseError("missing [ molecules ]");
    }
    return assemble_system(ir, molecule_types, system_spec);
}

SystemIR Parser::parse_from_file(const std::string& path) const
{
    std::vector<std::string> include_dirs;
    return parse_from_file(path, include_dirs);
}

SystemIR Parser::parse_from_file(const std::string& path,
                                 const std::vector<std::string>& include_dirs) const
{
    return parse_from_file(path, include_dirs, false);
}

SystemIR Parser::parse_from_file(const std::string& path,
                                 const std::vector<std::string>& include_dirs,
                                 bool dump_expanded) const
{
    try {
        std::string expanded = expand_topology(path, include_dirs, dump_expanded);
        return parse_from_string(expanded);
    } catch (const std::exception& e) {
        throw ParseError(e.what());
    }
}

} // namespace mdconv
