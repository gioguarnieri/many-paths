"""Catchment-sampled many-paths vulnerability index, per edge (v2).

Uniform OD sampling (trip_index.py) only defines the index on the 35-47%
of edges that sampled shortest paths happen to traverse -- the backbone
concentration of Kirkley et al. Here every edge gets its own OD sample
drawn from its *catchment*, mirroring Wang et al.'s finding that a road
segment's flow comes from a compact set of ~20 source zones.

Method, per directed edge e = (u -> v) of length l_e
----------------------------------------------------
1. Catchments: origins O = nodes that reach u (reverse Dijkstra), and
   destinations D = nodes reachable from v (forward Dijkstra), within
   MAX_TRIP_M - l_e of network distance. Interior nodes only.
2. Sample origins from O. One forward Dijkstra from each origin o gives
   the true shortest distance L0 to every candidate d in D at once.
   Accept (o, d) as a trip through e when
       MIN_TRIP_M <= T <= MAX_TRIP_M   and   T <= DELTA_REL * L0,
   where T = d(o,u) + l_e + d(v,d) is the through-e route length: the
   through route must sit inside the empirical tolerance envelope
   (Zhu & Levinson: most drivers accept routes within ~1.2-1.5x the
   shortest), not necessarily *be* the shortest -- that is what covers
   the low-betweenness loop edges uniform sampling misses.
   At most MAX_PER_ORIGIN destinations per origin, K_TRIPS trips total.
3. Suffering of each trip when e closes: greedy edge-disjoint o->d
   alternatives within budget DELTA_REL * L0, quality-weighted,
   s = exp(-lam * sum q_i)  (as trip_index.trip_suffering).
4. Each trip is weighted by its *reliance* on e,
       w = (DELTA_REL - T/L0) / (DELTA_REL - 1)  in (0, 1],
   1 when the through-e route is the trip's shortest path, decaying to 0
   at the tolerance limit -- the same linear-in-detour weighting used
   for alternative quality. Without it, marginal trips (whose shortest
   path survives the closure, capping their suffering at e^-lam) drown
   the dependent trips and compress the index: first unweighted runs
   gave 10th-90th percentile spans of only 0.35-0.56 (Sao Paulo) and
   0.29-0.35 (Tokyo).
5. V_e = weighted mean suffering over the edge's sampled trips; the CSV
   also reports the weight sum (effective dependent demand). Edges whose
   catchment yields no admissible trip are reported as no-demand.

Directed, giant strongly connected component, interior buffer: as v1.

Cost: all shortest-path work runs on a scipy.sparse.csgraph CSR core
(C Dijkstra, native cutoff via `limit`, no per-edge Python callbacks) --
benchmarked 4-15x faster than the tuned pure-Python A*/NetworkX paths and
the igraph/rustworkx alternatives on this workload. Edge closures poke the
collapsed CSR entry's weight (each entry keeps the multiset of its present
parallel-edge lengths, so closing one parallel exposes the next-cheapest,
exactly like removing that key from the MultiDiGraph). Trip sampling per
edge happens on the intact graph; closures apply only inside suffering
evaluation and are always restored.
"""

import argparse
import math
import random
import sys
import time
from bisect import insort
from pathlib import Path

import networkx as nx
import numpy as np
import osmnx as ox
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra as sp_dijkstra

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cities import CITIES, DEFAULT_RADIUS_M
from download import load_city
from redundancy import _haversine_m, BUFFER_M
from trip_index import DELTA_REL, MIN_TRIP_M, MAX_TRIP_M, LAM, MAX_ALTS

RESULTS_DIR = Path(__file__).resolve().parent.parent / "data" / "results"
FIG_DIR = Path(__file__).resolve().parent.parent / "figures"

K_TRIPS = 20         # target trips per edge (Wang et al.: ~20 source zones)
MAX_PER_ORIGIN = 2   # spread trips across origins, not one origin's fan-out
MAX_ORIGINS = 15     # origin sweeps per edge before giving up
SEED = 42


class CSRCore:
    """Directed street graph as a CSR matrix for scipy.sparse.csgraph.

    Parallel edges (u, v, *) collapse to one entry whose weight is the
    minimum length among the *currently open* parallels; `close`/`reopen`
    maintain that invariant, so pathfinding semantics match the
    MultiDiGraph with per-key removals. The reverse matrix is built once
    and used only on the intact graph (catchment sampling), never under
    closures.
    """

    def __init__(self, G):
        self.nodes = sorted(G.nodes)
        self.idx = {n: i for i, n in enumerate(self.nodes)}
        n = len(self.nodes)

        lengths = {}
        for u, v, d in G.edges(data=True):
            key = (self.idx[u], self.idx[v])
            lengths.setdefault(key, []).append(float(d.get("length", 1.0)))
        rows = np.fromiter((k[0] for k in lengths), dtype=np.int32,
                           count=len(lengths))
        cols = np.fromiter((k[1] for k in lengths), dtype=np.int32,
                           count=len(lengths))
        data = np.fromiter((min(v) for v in lengths.values()),
                           dtype=np.float64, count=len(lengths))
        self.A = csr_matrix((data, (rows, cols)), shape=(n, n))
        self.A_rev = self.A.transpose().tocsr()

        # slot of each collapsed edge in A.data + its open parallel lengths
        self.pos = {}
        indptr, indices = self.A.indptr, self.A.indices
        for ui in range(n):
            for j in range(indptr[ui], indptr[ui + 1]):
                self.pos[(ui, int(indices[j]))] = j
        self.avail = {}
        for key, ls in lengths.items():
            self.avail[self.pos[key]] = sorted(ls)

    def sweep(self, source, cutoff, reverse=False):
        """Distances (full float64 array, inf beyond cutoff)."""
        M = self.A_rev if reverse else self.A
        return sp_dijkstra(M, indices=source, limit=cutoff, directed=True)

    def query(self, source, target, budget):
        """Shortest path within budget as (dist, [node indices]), or None."""
        dist, pred = sp_dijkstra(self.A, indices=source,
                                 limit=np.nextafter(budget, np.inf),
                                 directed=True, return_predecessors=True)
        dv = dist[target]
        if not np.isfinite(dv) or dv > budget:
            return None
        path = [target]
        while path[-1] != source:
            path.append(int(pred[path[-1]]))
        path.reverse()
        return float(dv), path

    def close(self, ui, vi, length):
        """Close one parallel of (ui, vi); returns an undo token."""
        j = self.pos[(ui, vi)]
        self.avail[j].remove(length)
        self.A.data[j] = self.avail[j][0] if self.avail[j] else np.inf
        return j, length

    def close_min(self, ui, vi):
        """Close the currently cheapest parallel of (ui, vi)."""
        j = self.pos[(ui, vi)]
        return self.close(ui, vi, self.avail[j][0])

    def reopen(self, j, length):
        insort(self.avail[j], length)
        self.A.data[j] = self.avail[j][0]


def trip_suffering(core, oi, di, L0, closed, lam=LAM, max_alts=MAX_ALTS):
    """Suffering in [0,1] of trip (oi, di) when `closed` = (ui, vi, length)
    is removed. Same math as trip_index.trip_suffering."""
    budget = DELTA_REL * L0
    span = max(budget - L0, 1e-9)
    undo = [core.close(*closed)]
    quality_sum = 0.0
    n = 0
    while n < max_alts:
        alt = core.query(oi, di, budget)
        if alt is None:
            break
        dist, path = alt
        quality_sum += max(0.0, min(1.0, (budget - dist) / span))
        n += 1
        for a, b in zip(path[:-1], path[1:]):
            undo.append(core.close_min(a, b))
    for tok in reversed(undo):
        core.reopen(*tok)
    return math.exp(-lam * quality_sum)


def _ordered_within(dist, int_idx, radius):
    """Interior node indices with dist <= radius, ordered by distance
    (stable, so deterministic under ties). Returns a plain list."""
    di = dist[int_idx]
    mask = di <= radius
    sel = int_idx[mask]
    return sel[np.argsort(di[mask], kind="stable")].tolist()


def edge_catchment_trips(core, e, int_idx, rng, rev_cache):
    """Sample up to K_TRIPS trips (oi, di, L0, w) for interior edge
    e = (ui, vi, length), reliance weight w in (0, 1]."""
    ui, vi, l_e = e
    reach = max(MAX_TRIP_M - l_e, 0.0)

    # reverse sweep depends only on ui; edges arrive sorted by ui, so a
    # one-entry cache absorbs the ~2-3 out-edges sharing each tail
    if rev_cache.get("src") != ui:
        rev_cache["src"] = ui
        rev_cache["dist"] = core.sweep(ui, MAX_TRIP_M, reverse=True)
    dist_to_u = rev_cache["dist"]
    dist_from_v = core.sweep(vi, reach)

    origins = _ordered_within(dist_to_u, int_idx, reach)
    dests = _ordered_within(dist_from_v, int_idx, reach)
    if not origins or not dests:
        return []
    rng.shuffle(origins)

    trips = []
    for oi in origins[:MAX_ORIGINS]:
        d_ou = dist_to_u[oi]
        L0s = core.sweep(oi, MAX_TRIP_M)
        cands = rng.sample(dests, len(dests))
        n_from_o = 0
        for di in cands:
            if di == oi:
                continue
            T = d_ou + l_e + dist_from_v[di]
            if not (MIN_TRIP_M <= T <= MAX_TRIP_M):
                continue
            L0 = float(L0s[di])
            if not (0.0 < L0 < math.inf):
                continue
            w = (DELTA_REL - T / L0) / (DELTA_REL - 1.0)
            if w <= 0.0:
                continue
            trips.append((oi, di, L0, min(w, 1.0)))
            n_from_o += 1
            if n_from_o >= MAX_PER_ORIGIN or len(trips) >= K_TRIPS:
                break
        if len(trips) >= K_TRIPS:
            break
    return trips


def catchment_edge_vulnerability(G, interior_nodes, edges=None, seed=SEED,
                                 progress_every=500):
    """Returns ({edge: V_e}, {edge: (n_trips, weight_sum)}, [no-demand edges]).

    Edges are (u, v, key) in original node ids, as before.
    """
    rng = random.Random(seed)
    core = CSRCore(G)
    int_idx = np.fromiter(sorted(core.idx[n] for n in interior_nodes),
                          dtype=np.int64, count=len(interior_nodes))
    if edges is None:
        edges = [(u, v, k) for u, v, k in G.edges(keys=True)
                 if u in interior_nodes and v in interior_nodes]
    edges = sorted(edges)

    vuln, coverage, no_demand = {}, {}, []
    rev_cache = {}
    t0 = time.perf_counter()
    for i, (u, v, k) in enumerate(edges):
        ui, vi = core.idx[u], core.idx[v]
        l_e = float(G[u][v][k].get("length", 1.0))
        trips = edge_catchment_trips(core, (ui, vi, l_e), int_idx, rng,
                                     rev_cache)
        e = (u, v, k)
        if not trips:
            no_demand.append(e)
            continue
        w_sum = sum(w for _, _, _, w in trips)
        vuln[e] = sum(
            w * trip_suffering(core, oi, di, L0, (ui, vi, l_e))
            for oi, di, L0, w in trips) / w_sum
        coverage[e] = (len(trips), w_sum)
        if progress_every and (i + 1) % progress_every == 0:
            dt = time.perf_counter() - t0
            print(f"  ... {i + 1}/{len(edges)} edges "
                  f"[{dt:.0f}s, {1000 * dt / (i + 1):.0f} ms/edge]",
                  flush=True)
    return vuln, coverage, no_demand


def run_city(city, radius_m=DEFAULT_RADIUS_M, edge_limit=None, seed=SEED):
    import csv

    import matplotlib
    import matplotlib.pyplot as plt

    G = load_city(city, radius_m)
    giant = max(nx.strongly_connected_components(G), key=len)
    G = G.subgraph(giant).copy()
    clat, clon = CITIES[city]
    interior_nodes = {
        n for n, d in G.nodes(data=True)
        if _haversine_m(d["y"], d["x"], clat, clon) <= radius_m - BUFFER_M
    }
    interior_edges = [(u, v, k) for u, v, k in G.edges(keys=True)
                      if u in interior_nodes and v in interior_nodes]
    edges = interior_edges
    if edge_limit and edge_limit < len(edges):
        edges = random.Random(seed).sample(sorted(edges), edge_limit)

    t0 = time.perf_counter()
    vuln, coverage, no_demand = catchment_edge_vulnerability(
        G, interior_nodes, edges=edges, seed=seed)
    dt = time.perf_counter() - t0

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = RESULTS_DIR / f"{city}_{radius_m}m_catchment_vuln.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["u", "v", "key", "n_trips", "weight_sum", "vulnerability"])
        for (u, v, k), val in sorted(vuln.items()):
            n, ws = coverage[(u, v, k)]
            w.writerow([u, v, k, n, f"{ws:.2f}", f"{val:.4f}"])
        for u, v, k in sorted(no_demand):
            w.writerow([u, v, k, 0, "0.00", ""])

    cmap = matplotlib.colormaps["inferno"]
    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=1.0)
    colors, widths = [], []
    for u, v, k in G.edges(keys=True):
        e = (u, v, k)
        if e in vuln:
            rgba = cmap(norm(vuln[e]))
            colors.append(rgba)
            widths.append(0.7 + 2.0 * vuln[e])
        else:
            colors.append((0.35, 0.35, 0.35, 0.35))
            widths.append(0.4)
    fig, ax = ox.plot_graph(G, node_size=0, edge_color=colors,
                            edge_linewidth=widths, bgcolor="#111111",
                            show=False, close=False)
    sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=ax, shrink=0.7)
    cb.set_label("catchment vulnerability $V_e$", color="w")
    cb.ax.yaxis.set_tick_params(color="w")
    for lbl in cb.ax.get_yticklabels():
        lbl.set_color("w")
    ax.set_title(f"{city} catchment-sampled (directed, K={K_TRIPS}, "
                 f"$\\delta$={DELTA_REL})", color="w")
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / f"{city}_{radius_m}m_catchment_vuln.png", dpi=180,
                bbox_inches="tight", facecolor="#111111")
    plt.close(fig)

    n_scored = len(vuln)
    n_edges = len(edges)
    mean_v = sum(vuln.values()) / max(n_scored, 1)
    stranded = sum(1 for v_ in vuln.values() if v_ > 0.95)
    mean_trips = sum(n for n, _ in coverage.values()) / max(n_scored, 1)
    print(f"{city:12s} interior edges={n_edges:5,} scored={n_scored:5,} "
          f"({100 * n_scored / n_edges:.0f}%) no-demand={len(no_demand):4d} "
          f"mean trips/edge={mean_trips:.1f} mean V={mean_v:.3f} "
          f"stranding edges={stranded:4d} ({100 * stranded / max(n_scored, 1):.1f}%) "
          f"[{dt:.0f}s]")
    return vuln


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cities", nargs="*", default=None,
                    help="cities to run (default: all)")
    ap.add_argument("--limit", type=int, default=None,
                    help="score only a random subset of N interior edges")
    ap.add_argument("--radius", type=int, default=DEFAULT_RADIUS_M,
                    help="study radius in metres")
    args = ap.parse_args()
    for name in (args.cities or list(CITIES)):
        run_city(name, radius_m=args.radius, edge_limit=args.limit)
