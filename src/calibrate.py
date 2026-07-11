"""Calibration sweep for the v2 index parameters (gap #2).

Free parameters and their literature anchors:
- delta  (tolerance envelope): Xu & Chen elongation ratio 1.3-1.5 urban;
  Zhu & Levinson 70-80% of trips within +20% of shortest. Grid 1.1-1.6.
- lam    (EPD saturation): Rohrer & Sterbenz use lam=1. Grid 0.5/1/2.
- max_alts (portfolio cap): Zhu & Levinson observe ~2-3 routes per OD
  used interchangeably. Grid 2/3/4/6.

One pass covers the whole grid:
- greedy edge-disjoint extraction is nested -- the first A alternatives
  found with a_max=6 are exactly what a_max=A would find -- so per trip we
  record the quality prefix [q_1..q_m] once (per delta, since the budget
  delta*L0 shapes both which alts qualify and their q) and truncate;
- lam only enters post-hoc via V = sum(w * exp(-lam * sum q_i)) / sum(w);
- across deltas, trips share all Dijkstra sweeps and a common shuffled
  candidate order (common random numbers), so cross-delta differences are
  parameter effects, not sampling noise. Acceptance w > 0 is monotone in
  delta, so smaller-delta trip sets are near-subsets of larger-delta ones.

Outputs per city:
- data/results/calibration_{city}.csv: per (delta, lam, max_alts) summary
  (coverage, mean/std/deciles of V, stranding %, Spearman rank correlation
  vs the default combo delta=1.3/lam=1/A=4 on shared edges);
- data/results/calibration_{city}_rankstab.csv: delta x delta Spearman
  matrix at lam=1, A=4;
- figures/calibration_{city}.png: sensitivity panels.
"""

import argparse
import math
import random
import sys
import time
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cities import CITIES, DEFAULT_RADIUS_M
from download import load_city
from redundancy import _haversine_m, BUFFER_M
from trip_index import MIN_TRIP_M, MAX_TRIP_M
from catchment_index import (CSRCore, K_TRIPS, MAX_PER_ORIGIN, MAX_ORIGINS,
                             SEED, _ordered_within, RESULTS_DIR, FIG_DIR)

DELTAS = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
LAMS = [0.5, 1.0, 2.0]
ALTS = [2, 3, 4, 6]
A_MAX = max(ALTS)
BASE = (1.3, 1.0, 4)          # current defaults
N_EDGES = 1500


def edge_trips_multi(core, e, int_idx, rng, rev_cache):
    """Per-delta trip samples for one edge, sharing sweeps and candidate
    order. Returns {delta: [(oi, di, L0, w, T), ...]}."""
    ui, vi, l_e = e
    reach = max(MAX_TRIP_M - l_e, 0.0)

    if rev_cache.get("src") != ui:
        rev_cache["src"] = ui
        rev_cache["dist"] = core.sweep(ui, MAX_TRIP_M, reverse=True)
    dist_to_u = rev_cache["dist"]
    dist_from_v = core.sweep(vi, reach)

    origins = _ordered_within(dist_to_u, int_idx, reach)
    dests = _ordered_within(dist_from_v, int_idx, reach)
    trips = {dl: [] for dl in DELTAS}
    if not origins or not dests:
        return trips
    rng.shuffle(origins)

    for oi in origins[:MAX_ORIGINS]:
        if all(len(trips[dl]) >= K_TRIPS for dl in DELTAS):
            break
        d_ou = dist_to_u[oi]
        L0s = core.sweep(oi, MAX_TRIP_M)
        cands = rng.sample(dests, len(dests))
        n_from_o = {dl: 0 for dl in DELTAS}
        for di in cands:
            if di == oi:
                continue
            T = d_ou + l_e + dist_from_v[di]
            if not (MIN_TRIP_M <= T <= MAX_TRIP_M):
                continue
            L0 = float(L0s[di])
            if not (0.0 < L0 < math.inf):
                continue
            ratio = T / L0
            done = True
            for dl in DELTAS:
                w = (dl - ratio) / (dl - 1.0)
                if (w > 0.0 and n_from_o[dl] < MAX_PER_ORIGIN
                        and len(trips[dl]) < K_TRIPS):
                    trips[dl].append((oi, di, L0, min(w, 1.0), T))
                    n_from_o[dl] += 1
                if n_from_o[dl] < MAX_PER_ORIGIN and len(trips[dl]) < K_TRIPS:
                    done = False
            if done:
                break
    return trips


def trip_alt_qualities(core, oi, di, L0, closed, budget, a_max=A_MAX):
    """Quality prefix [q_1..q_m] of greedy edge-disjoint alternatives
    within budget, after closing `closed`. Truncating at A reproduces
    max_alts=A exactly (greedy nesting)."""
    span = max(budget - L0, 1e-9)
    undo = [core.close(*closed)]
    qs = []
    while len(qs) < a_max:
        alt = core.query(oi, di, budget)
        if alt is None:
            break
        dist, path = alt
        qs.append(max(0.0, min(1.0, (budget - dist) / span)))
        for a, b in zip(path[:-1], path[1:]):
            undo.append(core.close_min(a, b))
    for tok in reversed(undo):
        core.reopen(*tok)
    return qs


def calibrate_city(city, n_edges=N_EDGES, seed=SEED):
    import csv

    G = load_city(city, DEFAULT_RADIUS_M)
    G = G.subgraph(max(nx.strongly_connected_components(G), key=len)).copy()
    clat, clon = CITIES[city]
    interior_nodes = {
        n for n, d in G.nodes(data=True)
        if _haversine_m(d["y"], d["x"], clat, clon)
        <= DEFAULT_RADIUS_M - BUFFER_M
    }
    all_edges = sorted((u, v, k) for u, v, k in G.edges(keys=True)
                       if u in interior_nodes and v in interior_nodes)
    edges = random.Random(seed).sample(all_edges, min(n_edges, len(all_edges)))
    edges = sorted(edges)

    core = CSRCore(G)
    int_idx = np.fromiter(sorted(core.idx[n] for n in interior_nodes),
                          dtype=np.int64, count=len(interior_nodes))
    rng = random.Random(seed)
    rev_cache = {}

    # per delta: {edge: [(w, [q...]), ...]}
    samples = {dl: {} for dl in DELTAS}
    t0 = time.perf_counter()
    for i, (u, v, k) in enumerate(edges):
        ui, vi = core.idx[u], core.idx[v]
        l_e = float(G[u][v][k].get("length", 1.0))
        trips = edge_trips_multi(core, (ui, vi, l_e), int_idx, rng, rev_cache)
        for dl in DELTAS:
            recs = []
            for oi, di, L0, w, T in trips[dl]:
                qs = trip_alt_qualities(core, oi, di, L0, (ui, vi, l_e),
                                        budget=dl * L0)
                recs.append((w, qs))
            if recs:
                samples[dl][(u, v, k)] = recs
        if (i + 1) % 250 == 0:
            dt = time.perf_counter() - t0
            print(f"  ... {i + 1}/{len(edges)} edges "
                  f"[{dt:.0f}s, {1000 * dt / (i + 1):.0f} ms/edge]",
                  flush=True)

    def V_of(recs, lam, a):
        wsum = sum(w for w, _ in recs)
        return sum(w * math.exp(-lam * sum(qs[:a])) for w, qs in recs) / wsum

    def spearman(d1, d2):
        shared = sorted(set(d1) & set(d2))
        if len(shared) < 10:
            return float("nan")
        x = np.array([d1[e] for e in shared])
        y = np.array([d2[e] for e in shared])
        rx = np.argsort(np.argsort(x)).astype(float)
        ry = np.argsort(np.argsort(y)).astype(float)
        return float(np.corrcoef(rx, ry)[0, 1])

    # V maps per combo; base combo for rank comparisons
    Vmap = {}
    for dl in DELTAS:
        for lam in LAMS:
            for a in ALTS:
                Vmap[(dl, lam, a)] = {e: V_of(r, lam, a)
                                      for e, r in samples[dl].items()}
    base = Vmap[BASE]

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out = RESULTS_DIR / f"calibration_{city}.csv"
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["delta", "lam", "max_alts", "n_scored", "coverage_pct",
                    "mean_V", "std_V", "p10", "p50", "p90", "stranding_pct",
                    "spearman_vs_base", "mean_alts_found", "mean_sum_q"])
        for dl in DELTAS:
            nalts = [len(qs) for r in samples[dl].values() for _, qs in r]
            sumq = [sum(qs) for r in samples[dl].values() for _, qs in r]
            for lam in LAMS:
                for a in ALTS:
                    vm = Vmap[(dl, lam, a)]
                    vals = np.array(list(vm.values()))
                    w.writerow([
                        dl, lam, a, len(vals),
                        f"{100 * len(vals) / len(edges):.1f}",
                        f"{vals.mean():.4f}", f"{vals.std():.4f}",
                        f"{np.percentile(vals, 10):.4f}",
                        f"{np.percentile(vals, 50):.4f}",
                        f"{np.percentile(vals, 90):.4f}",
                        f"{100 * (vals > 0.95).mean():.2f}",
                        f"{spearman(vm, base):.4f}",
                        f"{np.mean(nalts):.2f}", f"{np.mean(sumq):.2f}",
                    ])

    stab = RESULTS_DIR / f"calibration_{city}_rankstab.csv"
    with open(stab, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["delta"] + [str(d) for d in DELTAS])
        for d1 in DELTAS:
            row = [spearman(Vmap[(d1, 1.0, 4)], Vmap[(d2, 1.0, 4)])
                   for d2 in DELTAS]
            w.writerow([d1] + [f"{r:.4f}" for r in row])

    _plot(city, samples, Vmap, edges)
    dt = time.perf_counter() - t0
    print(f"{city:12s} calibration on {len(edges):,} edges done [{dt:.0f}s] "
          f"-> {out.name}")


def _plot(city, samples, Vmap, edges):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    fig.suptitle(f"{city} — v2 parameter calibration "
                 f"({len(edges):,} edges, 5 km core)")

    ax = axes[0, 0]
    for lam, style in zip(LAMS, ["--", "-", ":"]):
        means = [np.mean(list(Vmap[(dl, lam, 4)].values())) for dl in DELTAS]
        ax.plot(DELTAS, means, style, marker="o", label=f"$\\lambda$={lam}")
    vals = [np.array(list(Vmap[(dl, 1.0, 4)].values())) for dl in DELTAS]
    ax.fill_between(DELTAS, [np.percentile(v, 10) for v in vals],
                    [np.percentile(v, 90) for v in vals], alpha=0.15,
                    label="p10–p90 ($\\lambda$=1)")
    ax.set_xlabel("$\\delta$"); ax.set_ylabel("V")
    ax.set_title("mean V vs tolerance (max_alts=4)")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    cov = [100 * len(samples[dl]) / len(edges) for dl in DELTAS]
    strand = [100 * (np.array(list(Vmap[(dl, 1.0, 4)].values())) > 0.95).mean()
              for dl in DELTAS]
    ax.plot(DELTAS, cov, marker="o", label="coverage %")
    ax.set_xlabel("$\\delta$"); ax.set_ylabel("coverage %")
    ax2 = ax.twinx()
    ax2.plot(DELTAS, strand, marker="s", color="tab:red", label="stranding %")
    ax2.set_ylabel("stranding %", color="tab:red")
    ax.set_title("coverage and stranding vs tolerance")

    ax = axes[1, 0]
    mat = np.array([[np.nan] * len(DELTAS) for _ in DELTAS])
    for i, d1 in enumerate(DELTAS):
        for j, d2 in enumerate(DELTAS):
            shared = set(Vmap[(d1, 1.0, 4)]) & set(Vmap[(d2, 1.0, 4)])
            x = np.array([Vmap[(d1, 1.0, 4)][e] for e in shared])
            y = np.array([Vmap[(d2, 1.0, 4)][e] for e in shared])
            rx = np.argsort(np.argsort(x)).astype(float)
            ry = np.argsort(np.argsort(y)).astype(float)
            mat[i, j] = np.corrcoef(rx, ry)[0, 1]
    im = ax.imshow(mat, vmin=0.5, vmax=1.0, cmap="viridis")
    ax.set_xticks(range(len(DELTAS)), [str(d) for d in DELTAS])
    ax.set_yticks(range(len(DELTAS)), [str(d) for d in DELTAS])
    for i in range(len(DELTAS)):
        for j in range(len(DELTAS)):
            ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center",
                    fontsize=7,
                    color="w" if mat[i, j] < 0.85 else "k")
    fig.colorbar(im, ax=ax, shrink=0.8)
    ax.set_title("rank stability across $\\delta$ (Spearman)")

    ax = axes[1, 1]
    for lam, style in zip(LAMS, ["--", "-", ":"]):
        means = [np.mean(list(Vmap[(1.3, lam, a)].values())) for a in ALTS]
        ax.plot(ALTS, means, style, marker="o", label=f"$\\lambda$={lam}")
    ax.set_xlabel("max_alts"); ax.set_ylabel("mean V")
    ax.set_title("portfolio saturation ($\\delta$=1.3)")
    ax.legend(fontsize=8)

    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / f"calibration_{city}.png", dpi=160,
                bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cities", nargs="*", default=None)
    ap.add_argument("--edges", type=int, default=N_EDGES)
    args = ap.parse_args()
    for name in (args.cities or ["sao_paulo", "tokyo"]):
        calibrate_city(name, n_edges=args.edges)
