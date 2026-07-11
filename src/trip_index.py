"""Trip-level many-paths vulnerability index, per edge (v1).

The local u->v detour framing (redundancy.py) over-punishes one-way systems:
closing a segment matters to the *trips through it*, not to its endpoints.

Method
------
1. Sample M random origin-destination node pairs inside the interior buffer,
   with shortest network distance in [MIN_TRIP_M, MAX_TRIP_M] (empirical:
   most urban trips are short).
2. Each trip's shortest path p0 (length L0) "uses" its edges. Tolerance
   budget per trip is relative: B = DELTA_REL * L0 (routes people accept
   are within ~1.2-1.5x the best route -- Zhu & Levinson 2015; Xu & Chen's
   elongation ratio).
3. For each edge e on p0: close e, then greedily extract up to MAX_ALTS
   edge-disjoint alternative o->d paths with length <= B. Alternative i of
   length L_i gets quality
       q_i = (B - L_i) / (B - L0)  in [0, 1]
   (1 = as good as the original trip, 0 = barely tolerable).
   Trip's suffering from losing e:  s = exp(-lam * sum_i q_i)
   (0 = plenty of good independent alternatives, 1 = trip is stranded).
4. Edge index V_e = mean suffering over sampled trips through e.
   Edges never traversed by a sampled trip are reported uncovered.

Directed throughout: one-way restrictions respected. Graph restricted to the
giant strongly connected component so o->d reachability is guaranteed
before closures.
"""

import math
import random
import sys
import time
from collections import defaultdict
from pathlib import Path

import networkx as nx
import osmnx as ox

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cities import CITIES, DEFAULT_RADIUS_M
from download import load_city
from redundancy import _cheapest_key, _shortest_alt

RESULTS_DIR = Path(__file__).resolve().parent.parent / "data" / "results"
FIG_DIR = Path(__file__).resolve().parent.parent / "figures"

N_TRIPS = 800        # sampled OD pairs per city
MIN_TRIP_M = 400.0   # shortest trips: below this, walking dominates anyway
MAX_TRIP_M = 5000.0  # trip-length literature: most urban trips < 5 km
DELTA_REL = 1.3      # calibrated: rankings robust 1.2-1.4; 1.1 is a
                     # different (near-shortest-only) regime -- calibrate.py
LAM = 1.0            # calibrated: rankings invariant for lam in [0.5, 2]
MAX_ALTS = 4         # calibrated: never binds (mean disjoint alts ~1.3-1.7)
MAX_TRIPS_PER_EDGE = 25
SEED = 42


def sample_trips(G, interior_nodes, n_trips, rng):
    """Random OD pairs with network distance in [MIN_TRIP_M, MAX_TRIP_M].

    Returns list of (o, d, L0, p0_edges) where p0_edges is the list of
    (u, v, key) edges of the shortest path.
    """
    nodes = list(interior_nodes)
    trips = []
    attempts = 0
    while len(trips) < n_trips and attempts < n_trips * 20:
        attempts += 1
        o = rng.choice(nodes)
        dists, paths = nx.single_source_dijkstra(G, o, cutoff=MAX_TRIP_M,
                                                 weight="length")
        cands = [n for n, dist in dists.items()
                 if dist >= MIN_TRIP_M and n in interior_nodes]
        if not cands:
            continue
        d = rng.choice(cands)
        p0 = paths[d]
        edges = [(a, b, _cheapest_key(G, a, b)) for a, b in zip(p0[:-1], p0[1:])]
        trips.append((o, d, dists[d], edges))
    return trips


def trip_suffering(G, o, d, L0, closed_edge, lam=LAM, max_alts=MAX_ALTS):
    """Suffering in [0,1] of trip (o,d) when closed_edge is removed."""
    u, v, k = closed_edge
    budget = DELTA_REL * L0
    span = max(budget - L0, 1e-9)

    removed = [(u, v, k, G[u][v][k])]
    G.remove_edge(u, v, k)

    quality_sum = 0.0
    n = 0
    while n < max_alts:
        alt = _shortest_alt(G, o, d, budget)
        if alt is None:
            break
        dist, path = alt
        quality_sum += max(0.0, min(1.0, (budget - dist) / span))
        n += 1
        for a, b in zip(path[:-1], path[1:]):
            kk = _cheapest_key(G, a, b)
            removed.append((a, b, kk, G[a][b][kk]))
            G.remove_edge(a, b, kk)

    for a, b, kk, dd in removed:
        G.add_edge(a, b, key=kk, **dd)
    return math.exp(-lam * quality_sum)


def trip_edge_vulnerability(G, interior_nodes, n_trips=N_TRIPS, seed=SEED):
    """Returns ({edge: V_e}, {edge: n_trips_through}, n_trips_sampled)."""
    rng = random.Random(seed)
    trips = sample_trips(G, interior_nodes, n_trips, rng)

    edge_trips = defaultdict(list)
    for o, d, L0, edges in trips:
        for e in edges:
            edge_trips[e].append((o, d, L0))

    vuln, coverage = {}, {}
    for e, tlist in edge_trips.items():
        if len(tlist) > MAX_TRIPS_PER_EDGE:
            tlist = rng.sample(tlist, MAX_TRIPS_PER_EDGE)
        s = [trip_suffering(G, o, d, L0, e) for o, d, L0 in tlist]
        vuln[e] = sum(s) / len(s)
        coverage[e] = len(tlist)
    return vuln, coverage, len(trips)


def run_city(city, radius_m=DEFAULT_RADIUS_M, n_trips=N_TRIPS):
    import csv

    import matplotlib
    import matplotlib.pyplot as plt

    G = load_city(city, radius_m)
    giant = max(nx.strongly_connected_components(G), key=len)
    G = G.subgraph(giant).copy()
    from redundancy import _haversine_m, BUFFER_M
    clat, clon = CITIES[city]
    interior_nodes = {
        n for n, d in G.nodes(data=True)
        if _haversine_m(d["y"], d["x"], clat, clon) <= radius_m - BUFFER_M
    }

    t0 = time.perf_counter()
    vuln, coverage, n_sampled = trip_edge_vulnerability(G, interior_nodes,
                                                        n_trips=n_trips)
    dt = time.perf_counter() - t0

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = RESULTS_DIR / f"{city}_{radius_m}m_trip_vuln.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["u", "v", "key", "n_trips", "vulnerability"])
        for (u, v, k), val in vuln.items():
            w.writerow([u, v, k, coverage[(u, v, k)], f"{val:.4f}"])

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
    cb.set_label("trip-level vulnerability $V_e$", color="w")
    cb.ax.yaxis.set_tick_params(color="w")
    for lbl in cb.ax.get_yticklabels():
        lbl.set_color("w")
    ax.set_title(f"{city} trip-level (directed, {n_sampled} trips, "
                 f"$\\delta$={DELTA_REL})", color="w")
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / f"{city}_{radius_m}m_trip_vuln.png", dpi=180,
                bbox_inches="tight", facecolor="#111111")
    plt.close(fig)

    covered = len(vuln)
    n_edges = G.number_of_edges()
    mean_v = sum(vuln.values()) / covered
    stranded = sum(1 for v_ in vuln.values() if v_ > 0.95)
    print(f"{city:12s} edges={n_edges:5,} covered={covered:5,} "
          f"({100 * covered / n_edges:.0f}%) mean V={mean_v:.3f} "
          f"stranding edges={stranded:4d} ({100 * stranded / covered:.1f}%) "
          f"[{dt:.0f}s]")
    return vuln


if __name__ == "__main__":
    only = sys.argv[1:] or list(CITIES)
    for name in only:
        run_city(name)
