#include "mdconv/gmx/expander.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace mdconv {

namespace {

// Minimal utilities. Keep it boring and deterministic.
static std::string ltrim(const std::string& s)
{
    size_t i = 0;
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) ++i;
    return s.substr(i);
}

static std::string parent_dir(const std::string& path)
{
    size_t last = path.find_last_of("/\\");
    if (last == std::string::npos) return ".";
    if (last == 0) return path.substr(0, 1);
    return path.substr(0, last);
}

static std::string join_path(const std::string& dir, const std::string& rel)
{
    if (dir.empty() || dir == ".") return rel;
    char last = dir[dir.size() - 1];
    if (last == '/' || last == '\\') return dir + rel;
    return dir + "/" + rel;
}

// v1: no filesystem canonicalization. We just normalize slashes a bit.
// (Important for reproducibility and to avoid platform-specific behavior.)
static std::string normalize_path(const std::string& p)
{
    std::string out = p;
    for (size_t i = 0; i < out.size(); ++i) {
        if (out[i] == '\\') out[i] = '/';
    }
    // remove redundant "./" segments (lightweight)
    // NOTE: we do NOT resolve ".." here; we forbid it earlier.
    while (out.rfind("./", 0) == 0) out.erase(0, 2);
    return out;
}

static bool is_absolute_path(const std::string& p)
{
    if (p.empty()) return false;
    if (p[0] == '/' || p[0] == '\\') return true;
    if (p.size() >= 2 && std::isalpha(static_cast<unsigned char>(p[0])) && p[1] == ':') return true;
    return false;
}

static bool file_readable(const std::string& p)
{
    std::ifstream in(p.c_str(), std::ios::in);
    return static_cast<bool>(in);
}

static bool has_allowed_ext(const std::string& p)
{
    // allowed: .top .itp .mdp .prm
    const char* exts[] = { "top", "itp", "mdp", "prm" };
    size_t dot = p.find_last_of('.');
    if (dot == std::string::npos) return false;
    if (dot + 1 >= p.size()) return false;

    // extension must be last component
    if (p.find('/', dot) != std::string::npos) return false;

    std::string ext = p.substr(dot + 1);
    for (size_t i = 0; i < 4; ++i) {
        if (ext == exts[i]) return true;
    }
    return false;
}

} // namespace

Expander::Frame::Frame() : line_no(0) {}

Expander::Expander(const std::string& root_path)
{
    push_file_or_throw(root_path);
}

void Expander::add_include_dir(const std::string& dir)
{
    if (!dir.empty()) {
        m_include_dirs.push_back(normalize_path(dir));
    }
}

bool Expander::parse_include_quotes_only(const std::string& line, std::string& target)
{
    std::string s = ltrim(line);
    if (s.empty() || s[0] != '#') return false;
    size_t i = 1;
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) ++i;
    if (s.compare(i, 7, "include") != 0) return false;
    i += 7;
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) ++i;

    // v1: only quotes
    if (i >= s.size() || s[i] != '"') {
        // if it looks like an include but not quotes, fail loudly later
        if (i < s.size() && s[i] == '<') {
            target.clear();
            throw std::runtime_error("v1 supports only #include \"...\" (angle includes are disabled)");
        }
        return false;
    }

    size_t start = i + 1;
    size_t end = s.find('"', start);
    if (end == std::string::npos || end == start) {
        throw std::runtime_error("invalid #include syntax (missing closing quote)");
    }

    target = s.substr(start, end - start);
    return true;
}

void Expander::validate_include_target_or_throw(const std::string& target)
{
    if (target.empty()) {
        throw std::runtime_error("invalid include path: empty");
    }
    std::string t = normalize_path(target);

    if (!has_allowed_ext(t)) {
        throw std::runtime_error("invalid include path: extension not allowed (allowed: .top .itp .mdp .prm)");
    }
}

void Expander::push_file_or_throw(const std::string& path)
{
    Frame f;
    f.path = normalize_path(path);
    f.dir  = parent_dir(f.path);
    f.line_no = 0;

    f.in.open(f.path.c_str(), std::ios::in);
    if (!f.in) {
        die_with_trace("cannot open file: " + f.path);
    }

    m_stack.push_back(Frame());
    m_stack.back().path = f.path;
    m_stack.back().dir = f.dir;
    m_stack.back().line_no = f.line_no;
    m_stack.back().in.swap(f.in);
}

bool Expander::path_in_stack(const std::string& canonical_path) const
{
    for (size_t i = 0; i < m_stack.size(); ++i) {
        if (m_stack[i].path == canonical_path) return true;
    }
    return false;
}

std::string Expander::resolve_include_or_throw(const Frame& includer,
                                              const std::string& target) const
{
    validate_include_target_or_throw(target);
    const std::string t = normalize_path(target);

    std::vector<std::string> tried;
    tried.reserve(1 + m_include_dirs.size());

    if (is_absolute_path(t)) {
        tried.push_back(t);
    } else {
    // 1) includer dir
    tried.push_back(normalize_path(join_path(includer.dir, t)));

    // 2) -I dirs
    for (size_t i = 0; i < m_include_dirs.size(); ++i) {
        tried.push_back(normalize_path(join_path(m_include_dirs[i], t)));
    }
    }

    for (size_t i = 0; i < tried.size(); ++i) {
        if (file_readable(tried[i])) {
            return tried[i];
        }
    }

    std::ostringstream oss;
    oss << "include not found: \"" << target << "\"\n";
    oss << "  search:\n";
    for (size_t i = 0; i < tried.size(); ++i) {
        oss << "    " << tried[i] << "\n";
    }
    die_with_trace(oss.str());
    return ""; // unreachable
}

[[noreturn]] void Expander::die_with_trace(const std::string& msg) const
{
    std::ostringstream oss;
    oss << "error: " << msg;
    if (!msg.empty() && msg[msg.size() - 1] != '\n') oss << "\n";
    oss << "include trace (most recent last):\n";
    for (size_t i = 0; i < m_stack.size(); ++i) {
        oss << "  " << m_stack[i].path << ":" << m_stack[i].line_no << "\n";
    }
    throw std::runtime_error(oss.str());
}

bool Expander::getline(std::string& out_line)
{
    while (!m_stack.empty()) {
        Frame& cur = m_stack.back();

        std::string line;
        if (!std::getline(cur.in, line)) {
            m_stack.pop_back();
            continue;
        }

        ++cur.line_no;
        if (!line.empty() && line.back() == '\r') line.pop_back();

        // handle includes
        try {
            std::string target;
            if (parse_include_quotes_only(line, target)) {
                std::string resolved = resolve_include_or_throw(cur, target);
                if (path_in_stack(resolved)) {
                    die_with_trace("recursive include detected: " + resolved);
                }
                push_file_or_throw(resolved);
                continue; // next line from included file
            }
        } catch (const std::runtime_error& e) {
            // rethrow with trace
            die_with_trace(e.what());
        }

        out_line = line;
        return true;
    }
    return false;
}

} // namespace mdconv
