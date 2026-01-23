#pragma once

#include <cstddef>
#include <fstream>
#include <string>
#include <vector>

namespace mdconv {

// Expands `#include "file.itp"` recursively and yields a flat stream of lines.
// v1 contract:
//   - only quotes includes: #include "path"
//   - allowed extensions: .top .itp .mdp .prm
//   - include target must be relative (no leading '/', no '..' component)
//   - search order: includer-dir first, then include_dirs (in order)
//   - fails loudly with include trace
class Expander {
public:
    explicit Expander(const std::string& root_path);

    // Add a directory to the -I search path.
    void add_include_dir(const std::string& dir);

    // Returns true and sets out_line if a line is available; false on EOF.
    // Errors throw std::runtime_error with an include trace.
    bool getline(std::string& out_line);

private:
    struct Frame {
        std::string   path;   // normalized resolved path
        std::string   dir;    // parent dir
        size_t        line_no;
        std::ifstream in;

        Frame();
    };

    std::vector<Frame> m_stack;
    std::vector<std::string> m_include_dirs;

private:
    static bool parse_include_quotes_only(const std::string& line, std::string& target);
    static void validate_include_target_or_throw(const std::string& target);

    void push_file_or_throw(const std::string& path);
    bool path_in_stack(const std::string& canonical_path) const;

    std::string resolve_include_or_throw(const Frame& includer,
                                         const std::string& target) const;

    [[noreturn]] void die_with_trace(const std::string& msg) const;
};

} // namespace mdconv
