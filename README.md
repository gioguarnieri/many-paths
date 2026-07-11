# many-paths

A pointwise (per-edge) vulnerability index for street networks, built on the
idea that a connection is less vulnerable when **many, mutually independent,
not-too-long paths** exist for the **trips that rely on it**.

## Catchment-sampled trip index (v2, current)

For each directed edge `e = (u → v)` of length `ℓ` (see
`src/catchment_index.py` for details):

1. **Catchment sampling**: origins that reach `u` and destinations reachable
   from `v` within network distance; accept an OD pair as a trip through `e`
   when the through-`e` route is inside the empirical tolerance envelope,
   `T = d(o,u) + ℓ + d(v,d) ≤ δ·L0` (δ = 1.3), with trip length in
   [0.4, 5] km. ~20 trips per edge, ≤2 per origin (Wang et al.'s ~20-zone
   catchment finding).
2. **Suffering per trip** when `e` closes: greedily extract edge-disjoint
   `o→d` alternatives within budget `δ·L0`, each weighted by quality
   `q_i = (budget − L_i)/(budget − L0) ∈ [0, 1]`; suffering
   `s = exp(−λ Σ q_i)` (EPD-style saturation, λ = 1).
3. **Reliance weighting**: each trip weighted by
   `w = (δ − T/L0)/(δ − 1) ∈ (0, 1]` — 1 when the through-`e` route is the
   shortest, 0 at the tolerance limit.
4. **Index**: `V_e` = reliance-weighted mean suffering. Coverage is ~100% of
   interior edges (uniform OD sampling only reaches 35–47%).

Directed throughout (one-way restrictions respected), giant strongly
connected component, 6 km download with 1 km boundary buffer (5 km scored
core). All shortest-path work runs on a `scipy.sparse.csgraph` CSR core;
temporary closures poke collapsed-entry weights with exact parallel-edge
semantics.

Grounding: stretch-bounded routes follow Xu & Chen (2018, elongation ratio
1.3–1.5 urban); saturating disjoint-path aggregation follows Rohrer &
Sterbenz's Effective Path Diversity; tolerance envelope, trip locality and
portfolio size from empirical mobility (Zhu & Levinson 2015; trip-length
literature). Full annotated bibliography and computational findings log:
`notes/LITERATURE.md`.

Earlier prototypes kept for comparison: `src/redundancy.py` (v0, local
endpoint detour — saturates in one-way systems) and `src/trip_index.py`
(v1, uniform OD sampling — coverage-limited).

## Layout

- `src/cities.py` — 10 most populous metro areas + central coordinates
- `src/download.py` — fetch & cache central drive networks (`python src/download.py [radius_m]`)
- `src/catchment_index.py` — the v2 index: CSV + map per city
- `src/redundancy.py`, `src/trip_index.py` — earlier prototypes (v0, v1)
- `notes/LITERATURE.md` — annotated bibliography + findings log
- `data/` (gitignored) — cached graphs and result CSVs
- `figures/` (gitignored) — vulnerability maps

## Setup

```
pip install -r requirements.txt
python src/download.py
python src/catchment_index.py sao_paulo tokyo
```
