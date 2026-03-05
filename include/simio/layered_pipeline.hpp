#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#include "simio/binning/FrameBinner.hpp"
#include "simio/intrinsics/WaterDipole.hpp"
#include "simio/regions/ChannelRegion.hpp"
#include "simio/simio.hpp"
namespace simio::layered {

enum class Region : uint8_t { Channel = 0, Reservoir = 1, Outside = 2 };

enum class Field : uint64_t {
    Layer0Bins = 1ull << 0,
    Layer0Region = 1ull << 1,
    KeyWrapped = 1ull << 2,
    KeyTimeWrapped = 1ull << 3,
    KeyCont = 1ull << 4,
    Velocity = 1ull << 5,
    Dipole = 1ull << 6,
    InChannel = 1ull << 7,
    TimeStepRaw = 1ull << 8,
    TimeStepMinImage = 1ull << 9,
    WrapK = 1ull << 10,
};

inline uint64_t field_mask(Field f) { return static_cast<uint64_t>(f); }

struct Layer0Config {
    int nx = 1;
    int nz = 1;

    // Optional wrapped centers; NaN means "box center".
    double x_center_wrapped = std::numeric_limits<double>::quiet_NaN();
    double z_center_wrapped = std::numeric_limits<double>::quiet_NaN();

    // Channel extents centered at (x_center_wrapped, z_center_wrapped).
    // A non-positive extent disables that axis constraint.
    double channel_length_x = 0.0;
    double channel_height_z = 0.0;

    // Gate planes expressed in center-relative z coordinates.
    // -1: z_rel < gate_left_rel_z, 0: between planes, +1: z_rel > gate_right_rel_z
    double gate_left_rel_z = 0.0;
    double gate_right_rel_z = 0.0;
};

struct Layer0Data {
    int nx = 1;
    int nz = 1;
    double Lx = 1.0;
    double Lz = 1.0;
    double dx = 1.0;
    double dz = 1.0;
    double x_center = 0.0;
    double z_center = 0.0;
    double channel_height_z = 0.0;

    std::vector<int> ix;
    std::vector<int> iz;
    std::vector<int> bin_id;
    std::vector<Region> region;
    // -1: left gate side, 0: between gates, +1: right gate side
    std::vector<int8_t> gate_side;

    void resize(size_t nmol) {
        ix.resize(nmol, 0);
        iz.resize(nmol, 0);
        bin_id.resize(nmol, 0);
        region.resize(nmol, Region::Outside);
        gate_side.resize(nmol, 0);
    }
};

class Layer0Geometry {
  public:
    explicit Layer0Geometry(const Layer0Config& cfg = {}) : cfg_(cfg) {}

    void set_config(const Layer0Config& cfg) { cfg_ = cfg; }
    const Layer0Config& config() const { return cfg_; }

    void compute(const Pbc3D& pbc, const std::vector<Vec3d>& key_wrapped, Layer0Data& out) const {
        binning::FrameBinnerConfig bcfg;
        bcfg.nx = cfg_.nx;
        bcfg.nz = cfg_.nz;
        bcfg.x_center_wrapped = cfg_.x_center_wrapped;
        bcfg.z_center_wrapped = cfg_.z_center_wrapped;

        regions::ChannelRegionConfig rcfg;
        rcfg.channel_length_x = cfg_.channel_length_x;
        rcfg.channel_height_z = cfg_.channel_height_z;
        rcfg.gate_left_rel_z = cfg_.gate_left_rel_z;
        rcfg.gate_right_rel_z = cfg_.gate_right_rel_z;

        binning::FrameBinner binner(bcfg);
        regions::ChannelRegion region_classifier(rcfg);
        const binning::FrameBinnerMeta meta = binner.make_meta(pbc);

        out.nx = meta.nx;
        out.nz = meta.nz;
        out.Lx = meta.Lx;
        out.Lz = meta.Lz;
        out.dx = meta.dx;
        out.dz = meta.dz;
        out.x_center = meta.x_center;
        out.z_center = meta.z_center;
        out.channel_height_z = std::max(0.0, cfg_.channel_height_z);

        out.resize(key_wrapped.size());

        for (size_t mid = 0; mid < key_wrapped.size(); ++mid) {
            const Vec3d& r = key_wrapped[mid];
            int ix = 0;
            int iz = 0;
            int bin = 0;
            binner.map_xz(pbc, meta, r, ix, iz, bin);

            const double x_rel = region_classifier.x_rel(pbc, r.v[0], out.x_center);
            const double z_rel = region_classifier.z_rel(pbc, r.v[2], out.z_center);

            out.ix[mid] = ix;
            out.iz[mid] = iz;
            out.bin_id[mid] = bin;
            out.region[mid] = region_classifier.in_channel(x_rel, z_rel) ? Region::Channel : Region::Reservoir;
            out.gate_side[mid] = region_classifier.gate_side(z_rel);
        }
    }

  private:
    Layer0Config cfg_{};
};

struct Layer1Config {
    // If true, water time tracking uses O position instead of COM.
    bool time_key_uses_ref_for_water = false;
};

struct Layer1Data {
    std::vector<Vec3d> key_wrapped;
    std::vector<Vec3d> key_time_wrapped;
    std::vector<Vec3d> ref_wrapped;
    std::vector<Vec3d> key_cont;
    std::vector<Vec3d> velocity;
    std::vector<Vec3d> dr_raw;
    std::vector<Vec3d> dr_mi;
    std::vector<std::array<int8_t, 3>> wrap_k;
    std::vector<Vec3d> dipole;
    std::vector<uint8_t> has_dipole;
    std::vector<uint8_t> is_in_channel;

    void resize(size_t nmol) {
        key_wrapped.resize(nmol, Vec3d{});
        key_time_wrapped.resize(nmol, Vec3d{});
        ref_wrapped.resize(nmol, Vec3d{});
        key_cont.resize(nmol, Vec3d{});
        velocity.resize(nmol, Vec3d{});
        dr_raw.resize(nmol, Vec3d{});
        dr_mi.resize(nmol, Vec3d{});
        wrap_k.resize(nmol, std::array<int8_t, 3>{0, 0, 0});
        dipole.resize(nmol, Vec3d{});
        has_dipole.resize(nmol, 0);
        is_in_channel.resize(nmol, 0);
    }
};

struct Layer1State {
    std::vector<uint8_t> has_prev;
    std::vector<Vec3d> prev_key_time_wrapped;
    std::vector<Vec3d> key_cont;

    bool has_prev_frame_time = false;
    double prev_frame_time_ps = 0.0;

    void resize(size_t nmol) {
        has_prev.resize(nmol, 0);
        prev_key_time_wrapped.resize(nmol, Vec3d{});
        key_cont.resize(nmol, Vec3d{});
    }
};

class Layer1Intrinsics {
  public:
    explicit Layer1Intrinsics(const Layer1Config& cfg = {}) : cfg_(cfg) {}

    void set_config(const Layer1Config& cfg) { cfg_ = cfg; }
    const Layer1Config& config() const { return cfg_; }

    void compute_local(const Topology& topo, const Frame& fr, Layer1Data& out, bool need_dipole) const {
        const size_t nmol = topo.mols.size();
        out.resize(nmol);

        for (size_t mid = 0; mid < nmol; ++mid) {
            const MolSpan& m = topo.mols[mid];

            Vec3d key{};
            Vec3d ref{};
            Vec3d mu{};
            bool has_ref = false;
            bool has_mu = false;

            if (m.type == MolType::Water) {
                assert(m.natoms >= 3);
                const int O = m.first + 0;
                const int H1 = m.first + 1;
                const int H2 = m.first + 2;

                if (need_dipole) {
                    const intrinsics::WaterDipoleResult wi = intrinsics::compute_water_dipole(fr, O, H1, H2);
                    key = wi.com_wrapped;
                    ref = wi.ref_wrapped;
                    mu = wi.dipole;
                    has_ref = true;
                    has_mu = wi.has_dipole;
                } else {
                    Vec3d rO = fr.atoms.pos(O);
                    const Vec3d rH1 = fr.atoms.pos(H1);
                    const Vec3d rH2 = fr.atoms.pos(H2);
                    fr.pbc.wrap_pos3(rO);

                    const Vec3d com_u = water_com_pbc_unwrapped(fr.pbc, rO, rH1, rH2);
                    Vec3d com_w = com_u;
                    fr.pbc.wrap_pos3(com_w);

                    key = com_w;
                    ref = rO;
                    has_ref = true;
                    has_mu = false;
                }
            } else {
                if (m.natoms >= 1) {
                    Vec3d r = fr.atoms.pos(m.first);
                    fr.pbc.wrap_pos3(r);
                    key = r;
                    ref = r;
                    has_ref = true;
                }
            }

            out.key_wrapped[mid] = key;
            out.ref_wrapped[mid] = ref;
            out.dipole[mid] = mu;
            out.has_dipole[mid] = has_mu ? 1 : 0;
            out.key_time_wrapped[mid] = (cfg_.time_key_uses_ref_for_water && m.type == MolType::Water && has_ref)
                                            ? ref
                                            : key;
        }
    }

    void update_time_and_flags(const Topology& topo, const Frame& fr, const Layer0Data& l0,
                               Layer1State& state, Layer1Data& out) const {
        const size_t nmol = topo.mols.size();
        state.resize(nmol);

        double dt_ps = 0.0;
        if (state.has_prev_frame_time) {
            dt_ps = fr.time_ps - state.prev_frame_time_ps;
            if (dt_ps <= 0.0) dt_ps = 0.0;
        }

        for (size_t mid = 0; mid < nmol; ++mid) {
            const Vec3d key_t = out.key_time_wrapped[mid];

            if (!state.has_prev[mid]) {
                state.has_prev[mid] = 1;
                state.prev_key_time_wrapped[mid] = key_t;
                state.key_cont[mid] = key_t;
                out.velocity[mid] = Vec3d{};
                out.dr_raw[mid] = Vec3d{};
                out.dr_mi[mid] = Vec3d{};
                out.wrap_k[mid] = std::array<int8_t, 3>{0, 0, 0};
            } else {
                const Vec3d dr_raw = key_t - state.prev_key_time_wrapped[mid];
                const Vec3d dr_mi = fr.pbc.min_image(dr_raw);
                std::array<int8_t, 3> kxyz{0, 0, 0};
                for (int a = 0; a < 3; ++a) {
                    if (fr.pbc.L[a] <= 0.0) {
                        kxyz[(size_t)a] = 0;
                        continue;
                    }
                    const double diff = dr_raw.v[a] - dr_mi.v[a];
                    const long long k_round = std::llround(diff / fr.pbc.L[a]);
                    if (k_round < -127) {
                        kxyz[(size_t)a] = -127;
                    } else if (k_round > 127) {
                        kxyz[(size_t)a] = 127;
                    } else {
                        kxyz[(size_t)a] = (int8_t)k_round;
                    }
                }

                state.prev_key_time_wrapped[mid] = key_t;
                state.key_cont[mid] = state.key_cont[mid] + dr_mi;
                out.velocity[mid] = (dt_ps > 0.0) ? (dr_mi / dt_ps) : Vec3d{};
                out.dr_raw[mid] = dr_raw;
                out.dr_mi[mid] = dr_mi;
                out.wrap_k[mid] = kxyz;
            }

            out.key_cont[mid] = state.key_cont[mid];
            out.is_in_channel[mid] = (mid < l0.region.size() && l0.region[mid] == Region::Channel) ? 1 : 0;
        }

        state.has_prev_frame_time = true;
        state.prev_frame_time_ps = fr.time_ps;
    }

  private:
    Layer1Config cfg_{};
};

struct FrameContext {
    const Topology& topo;
    const Frame& frame;
    const Layer0Data& layer0;
    const Layer1Data& layer1;
};

class PropertyKernel {
  public:
    virtual ~PropertyKernel() = default;
    virtual const char* name() const = 0;
    virtual uint64_t requires() const = 0;
    virtual void compute_frame(const FrameContext& ctx) = 0;
    virtual void finalize() {}
};

struct DipoleFrameSummary {
    int n_channel = 0;
    int n_reservoir = 0;
    double mean_mu_channel = 0.0;
    double mean_mu_reservoir = 0.0;
};

class DipoleProperty final : public PropertyKernel {
  public:
    const char* name() const override { return "DipoleProperty"; }

    uint64_t requires() const override {
        return field_mask(Field::Dipole) | field_mask(Field::Layer0Region);
    }

    void compute_frame(const FrameContext& ctx) override {
        DipoleFrameSummary s{};

        for (size_t mid = 0; mid < ctx.topo.mols.size(); ++mid) {
            if (!ctx.layer1.has_dipole[mid]) continue;
            const double mu = norm3(ctx.layer1.dipole[mid]);
            if (ctx.layer0.region[mid] == Region::Channel) {
                s.mean_mu_channel += mu;
                ++s.n_channel;
            } else if (ctx.layer0.region[mid] == Region::Reservoir) {
                s.mean_mu_reservoir += mu;
                ++s.n_reservoir;
            }
        }

        if (s.n_channel > 0) s.mean_mu_channel /= (double)s.n_channel;
        if (s.n_reservoir > 0) s.mean_mu_reservoir /= (double)s.n_reservoir;
        frames_.push_back(s);
    }

    void finalize() override {}

    const std::vector<DipoleFrameSummary>& frames() const { return frames_; }

  private:
    std::vector<DipoleFrameSummary> frames_;
};

enum class FluxSelection : uint8_t {
    All = 0,
    InChannelNow = 1,
    InChannelBoth = 2,
};

struct FluxGatingConfig {
    FluxSelection selection = FluxSelection::InChannelNow;
};

struct FluxFrameSummary {
    int64_t n_left_center = 0;
    int64_t n_right_center = 0;
    int64_t dn_center = 0;
    int64_t n_left_wrap = 0;
    int64_t n_right_wrap = 0;
    int64_t dn_wrap = 0;
};

class FluxGatingProperty final : public PropertyKernel {
  public:
    explicit FluxGatingProperty(const FluxGatingConfig& cfg = {}) : cfg_(cfg) {}

    const char* name() const override { return "FluxGatingProperty"; }

    uint64_t requires() const override {
        return field_mask(Field::Layer0Bins) | field_mask(Field::KeyTimeWrapped) |
               field_mask(Field::InChannel) | field_mask(Field::WrapK);
    }

    void compute_frame(const FrameContext& ctx) override {
        const size_t nmol = ctx.topo.mols.size();
        ensure_state(nmol);

        FluxFrameSummary fs{};

        for (size_t mid = 0; mid < nmol; ++mid) {
            const bool in_now = (ctx.layer1.is_in_channel[mid] != 0);
            const bool in_prev = (prev_in_channel_[mid] != 0);
            const bool include = include_molecule(in_now, in_prev);

            const double x_rel_now =
                ctx.frame.pbc.wrap_delta(0, ctx.layer1.key_time_wrapped[mid].v[0] - ctx.layer0.x_center);

            if (has_prev_x_rel_[mid] && include) {
                const double x_rel_prev = prev_x_rel_[mid];
                if (x_rel_prev <= 0.0 && x_rel_now > 0.0) {
                    ++fs.n_right_center;
                } else if (x_rel_prev >= 0.0 && x_rel_now < 0.0) {
                    ++fs.n_left_center;
                }
            }

            if (include) {
                const int kx = (int)ctx.layer1.wrap_k[mid][0];
                if (kx > 0) {
                    fs.n_right_wrap += kx;
                } else if (kx < 0) {
                    fs.n_left_wrap += -kx;
                }
            }

            prev_x_rel_[mid] = x_rel_now;
            has_prev_x_rel_[mid] = 1;
            prev_in_channel_[mid] = in_now ? 1 : 0;
        }

        fs.dn_center = fs.n_right_center - fs.n_left_center;
        fs.dn_wrap = fs.n_right_wrap - fs.n_left_wrap;

        cumulative_.n_left_center += fs.n_left_center;
        cumulative_.n_right_center += fs.n_right_center;
        cumulative_.dn_center = cumulative_.n_right_center - cumulative_.n_left_center;
        cumulative_.n_left_wrap += fs.n_left_wrap;
        cumulative_.n_right_wrap += fs.n_right_wrap;
        cumulative_.dn_wrap = cumulative_.n_right_wrap - cumulative_.n_left_wrap;

        frames_.push_back(fs);
    }

    void finalize() override {}

    const std::vector<FluxFrameSummary>& frames() const { return frames_; }
    const FluxFrameSummary& cumulative() const { return cumulative_; }

  private:
    bool include_molecule(bool in_now, bool in_prev) const {
        if (cfg_.selection == FluxSelection::All) return true;
        if (cfg_.selection == FluxSelection::InChannelBoth) return in_now && in_prev;
        return in_now;
    }

    void ensure_state(size_t nmol) {
        if (prev_x_rel_.size() == nmol) return;
        prev_x_rel_.assign(nmol, 0.0);
        has_prev_x_rel_.assign(nmol, 0);
        prev_in_channel_.assign(nmol, 0);
        frames_.clear();
        cumulative_ = FluxFrameSummary{};
    }

    FluxGatingConfig cfg_{};
    std::vector<double> prev_x_rel_;
    std::vector<uint8_t> has_prev_x_rel_;
    std::vector<uint8_t> prev_in_channel_;
    std::vector<FluxFrameSummary> frames_;
    FluxFrameSummary cumulative_{};
};

class Pipeline {
  public:
    Pipeline(const Layer0Config& l0_cfg = {}, const Layer1Config& l1_cfg = {})
        : layer0_(l0_cfg), layer1_(l1_cfg) {}

    void set_layer0_config(const Layer0Config& cfg) { layer0_.set_config(cfg); }
    void set_layer1_config(const Layer1Config& cfg) { layer1_.set_config(cfg); }

    void add_property(std::unique_ptr<PropertyKernel> kernel) {
        if (!kernel) return;
        required_mask_ |= kernel->requires();
        kernels_.emplace_back(std::move(kernel));
    }

    void process_frame(const Topology& topo, const Frame& frame) {
        const bool need_dipole = (required_mask_ & field_mask(Field::Dipole)) != 0ull;
        // Layer 1 (local intrinsics) computes per-molecule local fields.
        layer1_.compute_local(topo, frame, layer1_data_, need_dipole);

        // Layer 0 (pure geometry) maps positions to bins/regions.
        layer0_.compute(frame.pbc, layer1_data_.key_wrapped, layer0_data_);

        // Layer 1 time update uses Layer 0 region flags.
        layer1_.update_time_and_flags(topo, frame, layer0_data_, layer1_state_, layer1_data_);

        // Layer 2 property kernels consume reusable arrays.
        const FrameContext ctx{topo, frame, layer0_data_, layer1_data_};
        const uint64_t provided = field_mask(Field::Layer0Bins) | field_mask(Field::Layer0Region) |
                                  field_mask(Field::KeyWrapped) | field_mask(Field::KeyTimeWrapped) |
                                  field_mask(Field::KeyCont) | field_mask(Field::Velocity) |
                                  field_mask(Field::Dipole) | field_mask(Field::InChannel) |
                                  field_mask(Field::TimeStepRaw) | field_mask(Field::TimeStepMinImage) |
                                  field_mask(Field::WrapK);

        for (auto& k : kernels_) {
            if ((k->requires() & ~provided) != 0ull) {
                throw std::runtime_error(std::string("Unmet requirements in kernel: ") + k->name());
            }
            k->compute_frame(ctx);
        }
    }

    void finalize() {
        for (auto& k : kernels_) k->finalize();
    }

    const Layer0Data& layer0_data() const { return layer0_data_; }
    const Layer1Data& layer1_data() const { return layer1_data_; }

  private:
    Layer0Geometry layer0_;
    Layer1Intrinsics layer1_;

    Layer0Data layer0_data_;
    Layer1Data layer1_data_;
    Layer1State layer1_state_;

    std::vector<std::unique_ptr<PropertyKernel>> kernels_;
    uint64_t required_mask_ = 0ull;
};

}  // namespace simio::layered
