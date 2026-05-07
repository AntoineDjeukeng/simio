#include "simio/runtime/run_config.hpp"

#include <array>
#include <cctype>
#include <fstream>
#include <initializer_list>
#include <optional>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace simio::runtime {
namespace {

bool has_suffix(const std::string& s, const std::string& suffix) {
  if (suffix.size() > s.size()) return false;
  return s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string trim_copy(const std::string& s) {
  size_t b = 0, e = s.size();
  while (b < e && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
  while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
  return s.substr(b, e - b);
}

std::string upper_copy(std::string s) {
  for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

std::set<std::string> parse_name_set(const std::vector<std::string>& names) {
  std::set<std::string> out;
  for (const auto& n : names) {
    const std::string t = upper_copy(trim_copy(n));
    if (!t.empty()) out.insert(t);
  }
  return out;
}

std::vector<std::string> split_csv_names(const std::string& csv) {
  std::vector<std::string> out;
  std::string cur;
  auto flush = [&]() {
    const std::string t = trim_copy(cur);
    if (!t.empty()) out.push_back(t);
    cur.clear();
  };
  for (char c : csv) {
    if (c == ',' || std::isspace(static_cast<unsigned char>(c))) {
      flush();
    } else {
      cur.push_back(c);
    }
  }
  flush();
  return out;
}

std::string read_text_file(const std::string& path) {
  std::ifstream ifs(path);
  if (!ifs) throw std::runtime_error("Failed to open JSON config: " + path);
  std::ostringstream oss;
  oss << ifs.rdbuf();
  return oss.str();
}

std::optional<std::string> json_find_string(const std::string& text, const std::string& key) {
  const std::regex re("\"" + key + "\"\\s*:\\s*\"([^\"]*)\"");
  std::smatch m;
  if (std::regex_search(text, m, re)) return m[1].str();
  return std::nullopt;
}

std::optional<std::vector<std::string>> json_find_string_array(const std::string& text,
                                                               const std::string& key) {
  const std::regex re("\"" + key + "\"\\s*:\\s*\\[([^\\]]*)\\]");
  std::smatch m;
  if (!std::regex_search(text, m, re)) return std::nullopt;

  std::vector<std::string> out;
  const std::string body = m[1].str();
  const std::regex item_re("\"([^\"]*)\"");
  for (std::sregex_iterator it(body.begin(), body.end(), item_re), end; it != end; ++it) {
    out.push_back((*it)[1].str());
  }
  return out;
}

std::optional<std::string> json_find_number_token(const std::string& text, const std::string& key) {
  const std::regex re("\"" + key + "\"\\s*:\\s*([-+]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][-+]?\\d+)?)");
  std::smatch m;
  if (std::regex_search(text, m, re)) return m[1].str();
  return std::nullopt;
}

std::optional<std::string> json_find_string_any(const std::string& text,
                                                std::initializer_list<const char*> keys) {
  for (const char* k : keys) {
    if (auto v = json_find_string(text, k)) return v;
  }
  return std::nullopt;
}

std::optional<std::vector<std::string>> json_find_string_array_any(
    const std::string& text,
    std::initializer_list<const char*> keys) {
  for (const char* k : keys) {
    if (auto v = json_find_string_array(text, k)) return v;
  }
  return std::nullopt;
}

std::optional<std::string> json_find_number_any(const std::string& text,
                                                std::initializer_list<const char*> keys) {
  for (const char* k : keys) {
    if (auto v = json_find_number_token(text, k)) return v;
  }
  return std::nullopt;
}

std::optional<bool> json_find_bool_any(const std::string& text,
                                       std::initializer_list<const char*> keys) {
  for (const char* k : keys) {
    const std::regex re(std::string("\"") + k + "\"\\s*:\\s*(true|false)", std::regex::icase);
    std::smatch m;
    if (std::regex_search(text, m, re)) {
      return upper_copy(m[1].str()) == "TRUE";
    }
  }
  return std::nullopt;
}

std::optional<bool> json_find_bool_or_int_any(const std::string& text,
                                              std::initializer_list<const char*> keys) {
  if (auto b = json_find_bool_any(text, keys)) return b;
  if (auto n = json_find_number_any(text, keys)) {
    try {
      return std::stod(*n) != 0.0;
    } catch (...) {
      return std::nullopt;
    }
  }
  return std::nullopt;
}

int parse_int_text(const std::string& s, const char* name) {
  try {
    return std::stoi(s);
  } catch (const std::exception&) {
    throw std::runtime_error(std::string("Invalid integer for ") + name + ": " + s);
  }
}

double parse_double_text(const std::string& s, const char* name) {
  try {
    return std::stod(s);
  } catch (const std::exception&) {
    throw std::runtime_error(std::string("Invalid number for ") + name + ": " + s);
  }
}

std::array<double, 3> parse_array3_numbers(const std::string& s, const char* name) {
  std::array<double, 3> out{0.0, 0.0, 0.0};
  std::istringstream iss(s);
  std::string tok;
  for (int i = 0; i < 3; ++i) {
    if (!std::getline(iss, tok, ',')) {
      throw std::runtime_error(std::string("Expected 3 values for ") + name + ": " + s);
    }
    out[static_cast<size_t>(i)] = parse_double_text(trim_copy(tok), name);
  }
  return out;
}

struct MiniTopologyData {
  struct TypeEntry {
    std::string name;
    int molecules = 0;
    int atoms_per_molecule = 0; // 0 when unavailable in JSON
  };
  std::array<double, 3> channel_min{0.0, 0.0, 0.0};
  std::array<double, 3> channel_max{0.0, 0.0, 0.0};
  bool has_channel = false;
  std::vector<TypeEntry> types;
};

MiniTopologyData parse_mini_topology_json(const std::string& path) {
  const std::string text = read_text_file(path);
  MiniTopologyData out;

  {
    const std::regex ch_re("\"channel\"\\s*:\\s*\\{([\\s\\S]*?)\\}", std::regex::icase);
    std::smatch m;
    if (std::regex_search(text, m, ch_re)) {
      const std::string ch = m[1].str();
      const std::regex min_re("\"(?:min|min_nm)\"\\s*:\\s*\\[([^\\]]+)\\]", std::regex::icase);
      const std::regex max_re("\"(?:max|max_nm)\"\\s*:\\s*\\[([^\\]]+)\\]", std::regex::icase);
      std::smatch min_m, max_m;
      if (std::regex_search(ch, min_m, min_re) && std::regex_search(ch, max_m, max_re)) {
        out.channel_min = parse_array3_numbers(min_m[1].str(), "channel.min|min_nm");
        out.channel_max = parse_array3_numbers(max_m[1].str(), "channel.max|max_nm");
        out.has_channel = true;
      }
    }
  }

  const std::regex obj_re("\\{[^\\{\\}]*\\}");
  for (std::sregex_iterator it(text.begin(), text.end(), obj_re), end; it != end; ++it) {
    const std::string obj = it->str();
    const auto name_v = json_find_string(obj, "name");
    const auto mols_v = json_find_number_token(obj, "molecules");
    if (!name_v || !mols_v) continue;

    MiniTopologyData::TypeEntry e;
    e.name = upper_copy(trim_copy(*name_v));
    e.molecules = parse_int_text(*mols_v, "types[].molecules");
    if (auto apm_v = json_find_number_token(obj, "atoms_per_molecule")) {
      e.atoms_per_molecule = parse_int_text(*apm_v, "types[].atoms_per_molecule");
    }
    out.types.push_back(e);
  }

  return out;
}

void apply_topology_counts_from_mini_json(RunConfig& c) {
  if (c.topology_json.empty()) return;

  const MiniTopologyData topo = parse_mini_topology_json(c.topology_json);
  const auto water_set = parse_name_set(c.water_names);
  const auto na_set = parse_name_set(c.na_names);
  const auto cl_set = parse_name_set(c.cl_names);

  int nsol = 0;
  int nna = 0;
  int ncl = 0;
  bool has_sol = false;
  bool has_na = false;
  bool has_cl = false;
  std::vector<RunConfig::MolBlock> blocks;

  for (const auto& kv : topo.types) {
    const std::string& name = kv.name;
    const int count = kv.molecules;
    std::optional<simio::MolType> mapped;
    if (water_set.count(name) > 0) {
      nsol += count;
      has_sol = true;
      mapped = simio::MolType::Water;
    }
    if (na_set.count(name) > 0) {
      nna += count;
      has_na = true;
      mapped = simio::MolType::Cation;
    }
    if (cl_set.count(name) > 0) {
      ncl += count;
      has_cl = true;
      mapped = simio::MolType::Anion;
    }

    simio::MolType t = mapped.value_or(simio::MolType::Other);
    int natoms_per_mol = kv.atoms_per_molecule;
    if (t == simio::MolType::Water) {
      if (natoms_per_mol == 0) natoms_per_mol = 3;
      if (natoms_per_mol != 3) {
        throw std::runtime_error("topology_json water molecule '" + name +
                                 "' has atoms_per_molecule=" +
                                 std::to_string(natoms_per_mol) +
                                 ", but this driver requires 3 (OHH)");
      }
    } else if (t == simio::MolType::Cation || t == simio::MolType::Anion) {
      if (natoms_per_mol == 0) natoms_per_mol = 1;
      if (natoms_per_mol != 1) {
        throw std::runtime_error("topology_json ion molecule '" + name +
                                 "' has atoms_per_molecule=" +
                                 std::to_string(natoms_per_mol) +
                                 ", but this driver requires 1");
      }
    } else {
      if (natoms_per_mol <= 0) {
        throw std::runtime_error(
            "topology_json has unmapped molecule '" + name +
            "' without atoms_per_molecule; cannot preserve atom offset");
      }
    }

    if (count <= 0) continue;
    if (!blocks.empty() && blocks.back().type == t &&
        blocks.back().natoms_per_mol == natoms_per_mol) {
      blocks.back().nmol += count;
    } else {
      blocks.push_back(RunConfig::MolBlock{t, count, natoms_per_mol});
    }
  }

  if (has_sol) c.nsol = nsol;
  if (has_na) c.nna = nna;
  if (has_cl) c.ncl = ncl;
  if (!blocks.empty()) {
    c.mol_blocks = std::move(blocks);
    c.has_mol_blocks = true;
  }

  if (c.use_channel_bounds_from_topology) {
    if (!topo.has_channel) {
      throw std::runtime_error("topology_json missing channel.min/channel.max arrays: " + c.topology_json);
    }
    c.xmin = topo.channel_min[0];
    c.xmax = topo.channel_max[0];
    c.zmin = topo.channel_min[2];
    c.zmax = topo.channel_max[2];
  }
}

} // namespace

bool RunConfig::looks_like_json_path(const std::string& s) {
  return has_suffix(s, ".json");
}

RunConfig RunConfig::load(int argc, char** argv) {
  if (argc < 2) return from_cli(argc, argv);

  const std::string arg1 = argv[1];
  if (arg1 == "--config") {
    if (argc < 3) throw std::runtime_error("Missing JSON config path after --config.");
    return from_json_file(argv[2]);
  }
  if (argc == 2 && looks_like_json_path(arg1)) {
    return from_json_file(arg1);
  }

  return from_cli(argc, argv);
}

RunConfig RunConfig::from_cli(int argc, char** argv) {
  RunConfig c;
  if (argc < 2) return c;

  c.xtc_path = argv[1];
  if (argc > 2) c.max_frames = parse_int_text(argv[2], "max_frames");
  if (argc > 3) c.nthreads = parse_int_text(argv[3], "threads");
  if (argc > 4) c.grid_cell_nm = parse_double_text(argv[4], "grid_cell_nm");
  if (argc > 5) c.nsol = parse_int_text(argv[5], "nsol");
  if (argc > 6) c.nna = parse_int_text(argv[6], "nna");
  if (argc > 7) c.ncl = parse_int_text(argv[7], "ncl");
  if (argc > 8) c.xmin = parse_double_text(argv[8], "xmin");
  if (argc > 9) c.xmax = parse_double_text(argv[9], "xmax");
  if (argc > 10) c.zmin = parse_double_text(argv[10], "zmin");
  if (argc > 11) c.zmax = parse_double_text(argv[11], "zmax");
  if (argc > 12) c.nx = parse_int_text(argv[12], "nx");
  if (argc > 13) c.r_cw = parse_double_text(argv[13], "r_cw");
  if (argc > 14) c.r_aw = parse_double_text(argv[14], "r_aw");
  if (argc > 15) c.r_oo = parse_double_text(argv[15], "r_oo");
  if (argc > 16) c.gating_selection = argv[16];
  if (argc > 17) c.out_dir = argv[17];
  if (argc > 18) c.jump_keep_frames = parse_int_text(argv[18], "jump_keep_frames");
  if (argc > 19) c.frame_begin = parse_int_text(argv[19], "frame_begin");
  if (argc > 20) c.frame_end = parse_int_text(argv[20], "frame_end");
  if (argc > 21) c.r_nacl = parse_double_text(argv[21], "r_nacl");

  // Positional mode has no nz argument: keep behavior explicit.
  c.nz = c.nx;
  return c;
}

RunConfig RunConfig::from_json_file(const std::string& json_path) {
  const std::string text = read_text_file(json_path);
  RunConfig c;

  if (auto v = json_find_string(text, "xtc_path")) c.xtc_path = *v;
  const auto max_frames_v = json_find_number_any(text, {"max_frames"});
  const auto frame_begin_v = json_find_number_any(text, {"frame_begin", "begin_frame", "start_frame"});
  const auto frame_end_v = json_find_number_any(text, {"frame_end", "end_frame", "stop_frame"});
  const auto threads_v = json_find_number_any(text, {"threads", "nthreads"});
  const auto grid_cell_v = json_find_number_any(text, {"grid_cell_nm"});
  const auto nsol_v = json_find_number_any(text, {"nsol"});
  const auto nna_v = json_find_number_any(text, {"nna"});
  const auto ncl_v = json_find_number_any(text, {"ncl"});
  const auto xmin_v = json_find_number_any(text, {"xmin"});
  const auto xmax_v = json_find_number_any(text, {"xmax"});
  const auto zmin_v = json_find_number_any(text, {"zmin"});
  const auto zmax_v = json_find_number_any(text, {"zmax"});
  const auto nx_v = json_find_number_any(text, {"nx"});
  const auto nz_v = json_find_number_any(text, {"nz", "density_z_bins"});
  const auto rcw_v = json_find_number_any(text, {"r_cw"});
  const auto raw_v = json_find_number_any(text, {"r_aw"});
  const auto roo_v = json_find_number_any(text, {"r_oo"});
  const auto rnacl_v = json_find_number_any(text, {"r_nacl", "r_nacl_cluster", "nacl_cutoff_nm"});
  const auto gating_s_v = json_find_string_any(text, {"gating_sel", "gating_selection"});
  const auto out_dir_v = json_find_string_any(text, {"out_dir", "out_prefix"});
  const auto jump_keep_v = json_find_number_any(text, {"jump_keep_frames"});
  const auto bound_layer_v = json_find_number_any(text, {"bound_layer_nm", "bound_delta_nm", "adsorption_delta_nm"});

  if (auto v = json_find_string_any(text, {"topology_json", "mini_topology_json", "setup_json"})) {
    c.topology_json = *v;
  }
  if (auto v = json_find_bool_or_int_any(text, {"use_channel_bounds_from_topology", "use_topology_bounds"})) {
    c.use_channel_bounds_from_topology = *v;
  }

  if (auto v = json_find_string_array_any(text, {"water_names"})) {
    c.water_names = *v;
  } else if (auto s = json_find_string_any(text, {"water_names"})) {
    c.water_names = split_csv_names(*s);
  }
  if (auto v = json_find_string_array_any(text, {"na_names"})) {
    c.na_names = *v;
  } else if (auto s = json_find_string_any(text, {"na_names"})) {
    c.na_names = split_csv_names(*s);
  }
  if (auto v = json_find_string_array_any(text, {"cl_names"})) {
    c.cl_names = *v;
  } else if (auto s = json_find_string_any(text, {"cl_names"})) {
    c.cl_names = split_csv_names(*s);
  }

  if (max_frames_v) c.max_frames = parse_int_text(*max_frames_v, "max_frames");
  if (frame_begin_v) c.frame_begin = parse_int_text(*frame_begin_v, "frame_begin");
  if (frame_end_v) c.frame_end = parse_int_text(*frame_end_v, "frame_end");
  if (threads_v) c.nthreads = parse_int_text(*threads_v, "threads");
  if (grid_cell_v) c.grid_cell_nm = parse_double_text(*grid_cell_v, "grid_cell_nm");
  if (nsol_v) c.nsol = parse_int_text(*nsol_v, "nsol");
  if (nna_v) c.nna = parse_int_text(*nna_v, "nna");
  if (ncl_v) c.ncl = parse_int_text(*ncl_v, "ncl");
  if (xmin_v) c.xmin = parse_double_text(*xmin_v, "xmin");
  if (xmax_v) c.xmax = parse_double_text(*xmax_v, "xmax");
  if (zmin_v) c.zmin = parse_double_text(*zmin_v, "zmin");
  if (zmax_v) c.zmax = parse_double_text(*zmax_v, "zmax");
  if (nx_v) c.nx = parse_int_text(*nx_v, "nx");
  if (rcw_v) c.r_cw = parse_double_text(*rcw_v, "r_cw");
  if (raw_v) c.r_aw = parse_double_text(*raw_v, "r_aw");
  if (roo_v) c.r_oo = parse_double_text(*roo_v, "r_oo");
  if (rnacl_v) c.r_nacl = parse_double_text(*rnacl_v, "r_nacl");
  if (gating_s_v) c.gating_selection = *gating_s_v;
  if (out_dir_v) c.out_dir = *out_dir_v;
  if (jump_keep_v) c.jump_keep_frames = parse_int_text(*jump_keep_v, "jump_keep_frames");
  if (bound_layer_v) c.bound_layer_nm = parse_double_text(*bound_layer_v, "bound_layer_nm");

  if (nz_v) {
    c.nz = parse_int_text(*nz_v, "nz");
  } else {
    // Keep legacy behavior when nz isn't provided in JSON.
    c.nz = c.nx;
  }

  apply_topology_counts_from_mini_json(c);
  return c;
}

void RunConfig::validate_or_die() const {
  auto die = [&](const std::string& msg) { throw std::runtime_error("RunConfig invalid: " + msg); };

  if (xtc_path.empty()) die("xtc_path is empty");
  if (frame_begin < 0) die("frame_begin must be >= 0");
  if (frame_end != -1 && frame_end <= frame_begin) die("frame_end must be -1 or > frame_begin");
  if (max_frames <= 0) die("max_frames must be > 0");
  if (nthreads <= 0) die("nthreads must be > 0");
  if (grid_cell_nm <= 0.0) die("grid_cell_nm must be > 0");
  if (nsol < 0 || nna < 0 || ncl < 0) die("nsol/nna/ncl must be >= 0");
  if (nx <= 0 || nz <= 0) die("nx/nz must be > 0");
  if (r_cw <= 0.0 || r_aw <= 0.0 || r_oo <= 0.0 || r_nacl <= 0.0) {
    die("r_cw/r_aw/r_oo/r_nacl must be > 0");
  }
  if (jump_keep_frames <= 0) die("jump_keep_frames must be > 0");
  if (bound_layer_nm < 0.0) die("bound_layer_nm must be >= 0");
  if (out_dir.empty()) die("out_dir must be non-empty");
  if (has_mol_blocks) {
    for (size_t i = 0; i < mol_blocks.size(); ++i) {
      const auto& b = mol_blocks[i];
      if (b.nmol < 0) die("mol_blocks[" + std::to_string(i) + "].nmol must be >= 0");
      if (b.natoms_per_mol <= 0) {
        die("mol_blocks[" + std::to_string(i) + "].natoms_per_mol must be > 0");
      }
    }
  }
}

} // namespace simio::runtime
