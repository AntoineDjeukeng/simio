#pragma once

#include <string>
#include <vector>

#include "mdconv/ir/ir.hpp"

namespace mdconv {

struct GroFrame {
    std::vector<Vec3> x_nm;
    Box box_nm;
};

GroFrame parse_gro_from_string(const std::string& text);
GroFrame parse_gro_from_file(const std::string& path);

} // namespace mdconv
