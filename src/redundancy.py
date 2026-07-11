"""First prototype of the many-paths vulnerability index, per edge.

Idea being tested
-----------------
A directed edge e = (u -> v) is *not* vulnerable if, were it closed,
travellers could still get from u to v through many, mutually independent,
not-too-long alternative routes -- respecting one-way restrictions.

For each edge we greedily extract edge-disjoint alternative u->v directed
paths in G - e, accepting only paths within an absolute detour budget:

    len(alt) <= len(e) + delta_m          (delta_m ~ tolerable extra metres)

Each accepted alternative i contributes a quality weight

    q_i = 1 - extra_i / delta_m  in (0, 1]   (extra_i = len(alt_i) - len(e))

so a zero-cost parallel street counts fully and a barely-acceptable detour
counts almost nothing. Redundancy and vulnerability of the edge:

    k_e = sum_i q_i
    V_e = exp(-lam * k_e)   in (0, 1]    (1 = no alternative at all)

This is Rohrer & Sterbenz's EPD aggregation with full disjointness (greedy
removal), applied with the stretch-bounded route notion of Xu & Chen -- but
purely local to the edge, demand-free, and computable per edge.

Directionality: the graph is the raw OSM MultiDiGraph. A two-way street
contributes two directed edges scored independently; a one-way pair is NOT
an alternative for its opposite direction.

Boundary handling: edges whose endpoints lie within BUFFER_M of the study
radius are flagged exterior -- their alternatives may be cropped by the
download clip, so summary statistics use interior edges only and maps draw
exterior edges dimmed.

Remaining simplifications: disjointness is edge-based (shared intermediate
nodes allowed); the u->v detour is a local proxy -- OD-pair aggregation
within a travel radius is the planned next iteration.
"""

import math
import sys
import time
from heapq import heappop, heappush
from itertools import count
from pathlib import Path

import networkx as nx
import osmnx as ox

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cities import CITIES, DEFAULT_RADIUS_M
from download import load_city

RESULTS_DIR = Path(__file__).resolve().parent.parent / "data" / "results"
FIG_DIR = Path(__file__).resolve().parent.parent / "figures"

DELTA_M = 500.0   # tolerable extra detour in metres
LAM = 1.0         # EPD-style saturation constant
MAX_ALTS = 8      # stop looking after this many disjoint alternatives
BUFFER_M = 1000.0  # boundary buffer: edges this close to the clip are exterior


def _haversine_m(lat1, lon1, lat2, lon2):
    r = 6_371_000.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = p2 - p1
    dl = math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


def _ensure_xy(G):
    """Local metric (x, y) in metres per node, cached on the graph.

    Equirectangular projection about the mean latitude; distortion over a
    few-km study disc is < 1e-3 relative, absorbed by the heuristic's
    safety factor.
    """
    xy = G.graph.get("_xy")
    if xy is None:
        r = 6_371_000.0
        lat0 = math.radians(
            sum(d["y"] for _, d in G.nodes(data=True)) / len(G))
        kx = r * math.cos(lat0)
        kr = math.pi / 180.0
        xy = {n: (kx * kr * d["x"], r * kr * d["y"])
              for n, d in G.nodes(data=True)}
        G.graph["_xy"] = xy
    return xy


def _astar(G, source, target, budget, xy):
    """Budget-pruned A* over edge 'length': shortest directed source->target
    path no longer than budget, as (dist, path), else None.

    Straight-line distance is an admissible heuristic (edge length >= chord
    between endpoints; 0.999 covers projection distortion), so distances are
    exact. Nodes with g + h > budget cannot lie on a tolerable path and are
    pruned -- the search never leaves the budget ellipse around (o, d).
    """
    tx, ty = xy[target]
    hcache = {}

    def h(n):
        v = hcache.get(n)
        if v is None:
            x, y = xy[n]
            v = 0.999 * math.hypot(x - tx, y - ty)
            hcache[n] = v
        return v

    if h(source) > budget:
        return None
    c = count()
    fringe = [(h(source), next(c), source, 0.0, None)]
    enqueued = {}   # node -> best g pushed so far
    explored = {}   # node -> parent, set when settled
    while fringe:
        _, _, n, g, parent = heappop(fringe)
        if n in explored:
            continue
        explored[n] = parent
        if n == target:
            path = [n]
            while parent is not None:
                path.append(parent)
                parent = explored[parent]
            path.reverse()
            return g, path
        for nbr, keydata in G[n].items():
            if nbr in explored:
                continue
            ng = g + min(dd.get("length", 1.0) for dd in keydata.values())
            hn = h(nbr)
            if ng + hn > budget:
                continue
            prev = enqueued.get(nbr)
            if prev is not None and prev <= ng:
                continue
            enqueued[nbr] = ng
            heappush(fringe, (ng + hn, next(c), nbr, ng, n))
    return None


def _shortest_alt(H, u, v, budget):
    """Shortest directed u->v path in H no longer than budget, or None."""
    if u not in H or v not in H:
        return None
    return _astar(H, u, v, budget, _ensure_xy(H))


def _cheapest_key(H, a, b):
    """Key of the minimum-length parallel edge a->b (the one Dijkstra used)."""
    items = H.get_edge_data(a, b)
    return min(items, key=lambda k: items[k].get("length", 1.0))


def edge_vulnerability(G, delta_m=DELTA_M, lam=LAM, max_alts=MAX_ALTS):
    """Compute the index for every edge of a MultiDiGraph.

    Returns {(u, v, key): record} with vulnerability in (0, 1].
    """
    results = {}
    for u, v, key, data in list(G.edges(keys=True, data=True)):
        elen = float(data.get("length", 1.0))
        budget = elen + delta_m

        removed = [(u, v, key, data)]
        G.remove_edge(u, v, key)

        quality_sum = 0.0
        n_alts = 0
        first_extra = None
        while n_alts < max_alts:
            alt = _shortest_alt(G, u, v, budget)
            if alt is None:
                break
            dist, path = alt
            extra = max(dist - elen, 0.0)
            if first_extra is None:
                first_extra = extra
            quality_sum += 1.0 - extra / delta_m
            n_alts += 1
            for a, b in zip(path[:-1], path[1:]):
                k = _cheapest_key(G, a, b)
                removed.append((a, b, k, G[a][b][k]))
                G.remove_edge(a, b, k)

        for a, b, k, d in removed:
            G.add_edge(a, b, key=k, **d)

        results[(u, v, key)] = {
            "length": elen,
            "n_alts": n_alts,
            "first_extra_m": first_extra,
            "redundancy": quality_sum,
            "vulnerability": math.exp(-lam * quality_sum),
        }
    return results


def interior_mask(G, city, radius_m, buffer_m=BUFFER_M):
    """True for edges whose both endpoints are well inside the study radius."""
    clat, clon = CITIES[city]
    inside = {
        n: _haversine_m(d["y"], d["x"], clat, clon) <= radius_m - buffer_m
        for n, d in G.nodes(data=True)
    }
    return {
        (u, v, k): inside[u] and inside[v]
        for u, v, k in G.edges(keys=True)
    }


def run_city(city, radius_m=DEFAULT_RADIUS_M, delta_m=DELTA_M):
    """Load a cached city, compute the index, save CSV + map figure."""
    import csv

    G = load_city(city, radius_m)
    t0 = time.perf_counter()
    res = edge_vulnerability(G, delta_m=delta_m)
    dt = time.perf_counter() - t0
    interior = interior_mask(G, city, radius_m)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = RESULTS_DIR / f"{city}_{radius_m}m_edge_vuln.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["u", "v", "key", "length_m", "n_alts", "first_extra_m",
                    "redundancy", "vulnerability", "interior"])
        for (u, v, k), r in res.items():
            w.writerow([u, v, k, f"{r['length']:.1f}", r["n_alts"],
                        "" if r["first_extra_m"] is None else f"{r['first_extra_m']:.1f}",
                        f"{r['redundancy']:.4f}", f"{r['vulnerability']:.4f}",
                        int(interior[(u, v, k)])])

    # map coloured by vulnerability (dark = safe, bright = vulnerable);
    # exterior (buffer-zone) edges are dimmed
    import matplotlib
    import matplotlib.pyplot as plt
    cmap = matplotlib.colormaps["inferno"]
    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=1.0)
    colors, widths = [], []
    for u, v, k in G.edges(keys=True):
        x = res[(u, v, k)]["vulnerability"]
        rgba = cmap(norm(x))
        if not interior[(u, v, k)]:
            rgba = (rgba[0], rgba[1], rgba[2], 0.15)
        colors.append(rgba)
        widths.append(0.5 + 1.8 * x)
    fig, ax = ox.plot_graph(G, node_size=0, edge_color=colors,
                            edge_linewidth=widths, bgcolor="#111111",
                            show=False, close=False)
    sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=ax, shrink=0.7)
    cb.set_label("vulnerability  $V_e = e^{-k_e}$", color="w")
    cb.ax.yaxis.set_tick_params(color="w")
    for lbl in cb.ax.get_yticklabels():
        lbl.set_color("w")
    ax.set_title(f"{city}  (directed, r={radius_m} m, $\\Delta$={delta_m:.0f} m)",
                 color="w")
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig_path = FIG_DIR / f"{city}_{radius_m}m_vuln.png"
    fig.savefig(fig_path, dpi=180, bbox_inches="tight", facecolor="#111111")
    plt.close(fig)

    inner = [r for e, r in res.items() if interior[e]]
    n_in = len(inner)
    mean_v = sum(r["vulnerability"] for r in inner) / n_in
    no_alt = sum(1 for r in inner if r["n_alts"] == 0)
    print(f"{city:12s} nodes={len(G):5,} edges={G.number_of_edges():5,} "
          f"interior={n_in:5,} mean V={mean_v:.3f} "
          f"no-alt={no_alt:4d} ({100 * no_alt / n_in:.1f}%)  [{dt:.1f}s]")
    return res


if __name__ == "__main__":
    only = sys.argv[1:] or list(CITIES)
    for name in only:
        run_city(name)
