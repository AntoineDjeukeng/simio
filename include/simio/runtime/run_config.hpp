#pragma once

#include "simio/simio.hpp"

#include <string>
#include <vector>

namespace simio::runtime {

struct RunConfig {
  // input
  std::string xtc_path;

  // execution window
  int max_frames{100};
  int frame_begin{0};
  int frame_end{-1};

  // runtime
  int nthreads{4};

  // system sizes
  int nsol{5896};
  int nna{110};
  int ncl{110};

  // geometry / binning
  double xmin{7.11};
  double xmax{12.89};
  double zmin{0.901};
  double zmax{1.801};
  int nx{100};
  int nz{100}; // positional CLI doesn't provide nz; from_cli sets nz=nx.

  // model radii / cutoffs used by coord_x / etc.
  double r_cw{0.35};
  double r_aw{0.38};
  double r_oo{0.35};

  // cell list
  double grid_cell_nm{0.5};

  // analysis options
  std::string gating_selection{"all"};
  std::string out_dir{"."};
  int jump_keep_frames{50};
  double bound_layer_nm{0.30};

  // topology extras
  std::string topology_json;
  bool use_channel_bounds_from_topology{false};
  std::vector<std::string> water_names{"SOL", "WAT", "HOH"};
  std::vector<std::string> na_names{"NA", "NA+", "SOD"};
  std::vector<std::string> cl_names{"CL", "CL-", "CLA"};

  // Optional ordered molecule blocks inferred from topology_json types[] order.
  // This preserves atom offsets for unmapped species (e.g., walls).
  struct MolBlock {
    simio::MolType type{simio::MolType::Other};
    int nmol{0};
    int natoms_per_mol{0};
  };
  std::vector<MolBlock> mol_blocks;
  bool has_mol_blocks{false};

  // Loaders
  static RunConfig load(int argc, char** argv);             // auto CLI vs JSON
  static RunConfig from_cli(int argc, char** argv);         // positional CLI
  static RunConfig from_json_file(const std::string& path); // config.json
  void validate_or_die() const;

  // utilities
  static bool looks_like_json_path(const std::string& s);
};

} // namespace simio::runtime
