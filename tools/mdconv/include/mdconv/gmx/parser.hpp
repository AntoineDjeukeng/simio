#pragma once

#include <string>
#include "mdconv/ir/ir.hpp"

namespace mdconv {

class Parser {
public:
    Parser() = default;

    // v1: parse a topology-like string blob into IR.
    // current behavior: throws on empty/whitespace input.
    SystemIR parse_from_string(const std::string& input) const;

    // convenience for CLI; reads the whole file then delegates
    SystemIR parse_from_file(const std::string& path) const;
    SystemIR parse_from_file(const std::string& path,
                             const std::vector<std::string>& include_dirs) const;
    SystemIR parse_from_file(const std::string& path,
                             const std::vector<std::string>& include_dirs,
                             bool dump_expanded) const;

private:
    static bool is_all_whitespace(const std::string& s);
};

} // namespace mdconv
