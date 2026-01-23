#pragma once

#include <iosfwd>
#include <string>
#include "mdconv/ir/ir.hpp"

namespace mdconv {

class Writer {
public:
    Writer() = default;

    // v1: emit a stub LAMMPS header.
    void write_lammps_stub(std::ostream& os, const SystemIR& ir) const;

    // convenience for CLI
    void write_lammps_stub_file(const std::string& path, const SystemIR& ir) const;

    // v1: emit LAMMPS data format for atom_style full.
    void write_lammps_data(const SystemIR& ir, std::ostream& out) const;

    enum class KspaceMode { PPPM, None };

    struct Options {
        double cut_ang{10.0};
        KspaceMode kspace{KspaceMode::PPPM};
        double pppm_accuracy{1e-4};
    };

    // v1: emit a LAMMPS input snippet for coefficients and styles.
    void write_lammps_input_snippet(const SystemIR& ir,
                                    std::ostream& out,
                                    const Options& options) const;
};

} // namespace mdconv
