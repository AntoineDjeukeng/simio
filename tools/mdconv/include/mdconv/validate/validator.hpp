#pragma once

#include "mdconv/ir/ir.hpp"

namespace mdconv {

class Validator {
public:
    Validator() = default;

    // v1: fails if zero atoms.
    void validate_or_throw(const SystemIR& ir) const;
};

} // namespace mdconv
