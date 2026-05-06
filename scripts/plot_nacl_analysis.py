#!/usr/bin/env python3
import argparse
import csv
import os


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def col(rows, name):
    return [float(r[name]) for r in rows]


def maybe_col(rows, name):
    if not rows or name not in rows[0]:
        return None
    return col(rows, name)


def add_channel_edges(ax, x_left, x_right):
    ax.axvline(x_left, color="0.25", lw=1, ls="--")
    ax.axvline(x_right, color="0.25", lw=1, ls="--")
    ax.axvspan(x_left, x_right, color="0.9", zorder=-10)


def savefig(fig, out):
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--x-left", type=float, default=8.048)
    ap.add_argument("--x-right", type=float, default=13.828)
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = os.path.join(args.results_dir, "figures")
    os.makedirs(out_dir, exist_ok=True)

    density = read_csv(os.path.join(args.results_dir, "density_x.csv"))
    assoc = read_csv(os.path.join(args.results_dir, "nacl_association_x.csv"))
    gating = read_csv(os.path.join(args.results_dir, "gating_flux.csv"))
    clusters = read_csv(os.path.join(args.results_dir, "nacl_cluster_x.csv"))

    x = col(density, "x_center_nm")

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(x, col(density, "rho_water_mean"), label="SOL")
    ax.plot(x, col(density, "rho_na_mean"), label="NA")
    ax.plot(x, col(density, "rho_cl_mean"), label="CL")
    add_channel_edges(ax, args.x_left, args.x_right)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("density")
    ax.set_title("Species density vs x")
    ax.legend()
    savefig(fig, os.path.join(out_dir, "density_x.png"))

    xa = col(assoc, "x_center_nm")
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(xa, col(assoc, "N_CIP_mean"), label="CIP")
    ax.plot(xa, col(assoc, "N_SSIP_mean"), label="SSIP")
    ax.plot(xa, col(assoc, "N_ASSOC_mean"), label="ASSOC")
    add_channel_edges(ax, args.x_left, args.x_right)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("pair count / frame")
    ax.set_title("Na-Cl association vs x")
    ax.legend()
    savefig(fig, os.path.join(out_dir, "nacl_association_counts_x.png"))

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(xa, col(assoc, "f_CIP_mean"), label="f_CIP")
    ax.plot(xa, col(assoc, "f_SSIP_mean"), label="f_SSIP")
    ax.plot(xa, col(assoc, "f_bridge_mean"), label="f_bridge")
    add_channel_edges(ax, args.x_left, args.x_right)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("fraction")
    ax.set_ylim(-0.05, 1.05)
    ax.set_title("Association fractions vs x")
    ax.legend()
    savefig(fig, os.path.join(out_dir, "nacl_association_fractions_x.png"))

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(xa, col(assoc, "CN_NaOW_mean"), label="CN_NaOW")
    ax.plot(xa, col(assoc, "CN_ClOW_mean"), label="CN_ClOW")
    ax.plot(xa, col(assoc, "N_bridge_water_mean"), label="bridge water")
    ax.plot(xa, col(assoc, "N_bridged_pair_mean"), label="bridged pair")
    add_channel_edges(ax, args.x_left, args.x_right)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("count / coordination")
    ax.set_title("Hydration and bridge metrics vs x")
    ax.legend()
    savefig(fig, os.path.join(out_dir, "hydration_bridge_x.png"))

    xc = col(clusters, "x_center_nm")
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(xc, col(clusters, "nacl_cluster_count_mean"), label="cluster count")
    ax.plot(xc, col(clusters, "nacl_cluster_size_mean"), label="cluster size")
    add_channel_edges(ax, args.x_left, args.x_right)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("mean")
    ax.set_title("Simple Na-Cl cluster graph vs x")
    ax.legend()
    savefig(fig, os.path.join(out_dir, "nacl_clusters_x.png"))

    tg = [float(r["time_ps"]) / 1000.0 for r in gating]
    fig, ax = plt.subplots(figsize=(8, 4))
    for sp in ("water", "na", "cl"):
        right = col(gating, f"center_{sp}_cum_right")
        left_events = col(gating, f"center_{sp}_left")
        left_cum = []
        s = 0.0
        for v in left_events:
            s += v
            left_cum.append(s)
        net = [r - l for r, l in zip(right, left_cum)]
        ax.plot(tg, net, label=sp.upper())
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("center net crossings")
    ax.set_title("Net center crossings")
    ax.legend()
    savefig(fig, os.path.join(out_dir, "gating_net_crossings.png"))

    print("sanity:")
    print(f"  frames in gating_flux: {len(gating)}")
    print(f"  x bins: {len(assoc)}")
    print(f"  final time ns: {tg[-1] if tg else 0.0}")
    print(f"  final water cum right: {gating[-1]['center_water_cum_right'] if gating else 'NA'}")
    print(f"  final na cum right: {gating[-1]['center_na_cum_right'] if gating else 'NA'}")
    print(f"  final cl cum right: {gating[-1]['center_cl_cum_right'] if gating else 'NA'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
