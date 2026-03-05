#pragma once

#include <cstdint>

namespace simio::properties {

enum class MoleculeSelection : uint8_t {
    All = 0,
    InChannelNow = 1,
    InChannelBoth = 2,
};

inline bool include_molecule(MoleculeSelection selection, bool in_now, bool in_prev) {
    if (selection == MoleculeSelection::All) return true;
    if (selection == MoleculeSelection::InChannelBoth) return in_now && in_prev;
    return in_now;
}

struct GatingFrameCount {
    int64_t n_left = 0;
    int64_t n_right = 0;
    int64_t dn = 0;
};

}  // namespace simio::properties

