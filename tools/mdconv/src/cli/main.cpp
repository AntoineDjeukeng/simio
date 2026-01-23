#include "mdconv/gmx/parser.hpp"
#include "mdconv/validate/validator.hpp"
#include "mdconv/backend/writer.hpp"
#include "mdconv/diag/errors.hpp"
#include "mdconv/util/dump.hpp"
#include "mdconv/gmx/gro.hpp"
#include "mdconv/normalize/resolver.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>

static void usage(const char* prog)
{
    std::cerr << "Usage:\n"
              << "  " << prog << " <input_top_like.txt> <output.data> [options]\n"
              << "Options:\n"
              << "  --gro <file.gro>\n"
              << "  -I <dir> (add include directory; can be repeated)\n"
              << "  --include <dir>\n"
              << "  --dump-expanded\n"
              << "  --cut <angstrom>\n"
              << "  --kspace pppm|none\n"
              << "  --pppm <accuracy>\n";
}

int main(int argc, char** argv)
{
    if (argc < 3) {
        usage(argv[0]);
        return 1;
    }

    const std::string input_path = argv[1];
    const std::string output_path = argv[2];
    std::string gro_path;
    std::vector<std::string> include_dirs;
    mdconv::Writer::Options wopts;
    bool pppm_set = false;
    bool dump_expanded = false;
    for (int i = 3; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--gro" && i + 1 < argc) {
            gro_path = argv[++i];
        } else if (arg == "-I" && i + 1 < argc) {
            include_dirs.push_back(argv[++i]);
        } else if (arg == "--include" && i + 1 < argc) {
            include_dirs.push_back(argv[++i]);
        } else if (arg.size() > 2 && arg.rfind("-I", 0) == 0) {
            include_dirs.push_back(arg.substr(2));
        } else if (arg == "--dump-expanded") {
            dump_expanded = true;
        } else if (arg == "--cut" && i + 1 < argc) {
            wopts.cut_ang = std::atof(argv[++i]);
        } else if (arg == "--kspace" && i + 1 < argc) {
            std::string mode = argv[++i];
            if (mode == "pppm") {
                wopts.kspace = mdconv::Writer::KspaceMode::PPPM;
            } else if (mode == "none") {
                wopts.kspace = mdconv::Writer::KspaceMode::None;
            } else {
                usage(argv[0]);
                return 1;
            }
        } else if (arg == "--pppm" && i + 1 < argc) {
            wopts.pppm_accuracy = std::atof(argv[++i]);
            pppm_set = true;
        } else {
            usage(argv[0]);
            return 1;
        }
    }
    if (wopts.cut_ang <= 0.0) {
        std::cerr << "error: --cut must be positive\n";
        return 1;
    }
    if (wopts.pppm_accuracy <= 0.0) {
        std::cerr << "error: --pppm must be positive\n";
        return 1;
    }
    if (wopts.kspace == mdconv::Writer::KspaceMode::None) {
        if (pppm_set) {
            std::cerr << "error: --pppm requires --kspace pppm\n";
            return 1;
        }
        std::cerr << "warning: kspace disabled; electrostatics will be truncated\n";
    }

    try {
        mdconv::Parser parser;
        mdconv::Validator validator;
        mdconv::Writer writer;

        mdconv::SystemIR ir = parser.parse_from_file(input_path, include_dirs, dump_expanded);
        if (!gro_path.empty()) {
            mdconv::GroFrame gro = mdconv::parse_gro_from_file(gro_path);
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
        ir = mdconv::resolve_system(ir);
        validator.validate_or_throw(ir);
        const char* dump_env = std::getenv("MDCONV_DUMP");
        if (dump_env != nullptr && std::string(dump_env) == "1") {
            std::cout << mdconv::to_json(ir) << "\n";
        }
        std::ofstream out(output_path.c_str(), std::ios::out | std::ios::trunc);
        if (!out) {
            throw std::runtime_error("failed to open output file: " + output_path);
        }
        writer.write_lammps_data(ir, out);
        writer.write_lammps_input_snippet(ir, std::cout, wopts);

        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 2;
    }
}
