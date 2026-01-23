#include "mdconv/util/dump.hpp"

#include <iomanip>
#include <sstream>

namespace mdconv {

static std::string json_escape(const std::string& s)
{
    std::ostringstream out;
    for (size_t i = 0; i < s.size(); ++i) {
        unsigned char c = static_cast<unsigned char>(s[i]);
        switch (c) {
        case '\\':
            out << "\\\\";
            break;
        case '"':
            out << "\\\"";
            break;
        case '\b':
            out << "\\b";
            break;
        case '\f':
            out << "\\f";
            break;
        case '\n':
            out << "\\n";
            break;
        case '\r':
            out << "\\r";
            break;
        case '\t':
            out << "\\t";
            break;
        default:
            if (c < 0x20) {
                out << "\\u"
                    << std::hex << std::setw(4) << std::setfill('0')
                    << static_cast<int>(c)
                    << std::dec << std::setw(0);
            } else {
                out << s[i];
            }
            break;
        }
    }
    return out.str();
}

static const char* unit_system_name(UnitSystem units)
{
    switch (units) {
    case UnitSystem::GMX:
        return "GMX";
    }
    return "unknown";
}

static const char* mixing_rule_name(MixingRule rule)
{
    switch (rule) {
    case MixingRule::LorentzBerthelot:
        return "LorentzBerthelot";
    }
    return "unknown";
}

std::string to_json(const SystemIR& ir)
{
    std::ostringstream out;
    out << std::setprecision(17);

    out << "{";
    out << "\"unit_system\":\"" << unit_system_name(ir.units) << "\",";

    out << "\"counts\":{";
    out << "\"atom_types\":" << ir.atom_types.size() << ",";
    out << "\"atoms\":" << ir.atoms.size();
    out << "},";

    out << "\"atom_types\":[";
    for (size_t i = 0; i < ir.atom_types.size(); ++i) {
        const AtomType& at = ir.atom_types[i];
        if (i != 0) {
            out << ",";
        }
        out << "{";
        out << "\"id\":" << at.id << ",";
        out << "\"name\":\"" << json_escape(at.name) << "\",";
        out << "\"mass\":" << at.mass_amu << ",";
        out << "\"sigma\":" << at.sigma_nm << ",";
        out << "\"epsilon\":" << at.epsilon_kj;
        out << "}";
    }
    out << "],";

    out << "\"atoms\":[";
    for (size_t i = 0; i < ir.atoms.size(); ++i) {
        const Atom& a = ir.atoms[i];
        if (i != 0) {
            out << ",";
        }
        out << "{";
        out << "\"id\":" << a.id << ",";
        out << "\"molecule_id\":" << a.molecule_id << ",";
        out << "\"type_id\":" << a.type_id << ",";
        out << "\"charge\":" << a.charge_e << ",";
        out << "\"x\":" << a.x_nm << ",";
        out << "\"y\":" << a.y_nm << ",";
        out << "\"z\":" << a.z_nm;
        out << "}";
    }
    out << "],";

    out << "\"nb_defaults\":{";
    out << "\"mixing\":\"" << mixing_rule_name(ir.nb_defaults.mixing) << "\",";
    out << "\"fudge_lj\":" << ir.nb_defaults.fudge_lj << ",";
    out << "\"fudge_qq\":" << ir.nb_defaults.fudge_qq;
    out << "},";

    out << "\"box\":{";
    out << "\"lx\":" << ir.box.lx_nm << ",";
    out << "\"ly\":" << ir.box.ly_nm << ",";
    out << "\"lz\":" << ir.box.lz_nm;
    out << "}";

    out << "}";
    return out.str();
}

} // namespace mdconv
