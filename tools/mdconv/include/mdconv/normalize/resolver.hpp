#pragma once

#include "mdconv/ir/ir.hpp"

namespace mdconv {

// v1: pass-through placeholder for semantic resolution.
SystemIR resolve_system(const SystemIR& ir);

} // namespace mdconv
