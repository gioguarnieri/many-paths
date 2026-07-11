# Literature notes — many-paths vulnerability index

Working notes with citations. Each claim carries: paper link + where in the
paper the info was found. Items marked ⚠ were taken from abstracts/search
snippets and still need full-text verification.

## 1. Closest competitor measures

### Xu & Chen (2018) — Transportation network redundancy
- Journal: [Transportation Research Part B 114:68–85](https://www.sciencedirect.com/science/article/abs/pii/S0191261517306781);
  full-text open report version: [MPC-17-327, UGPTI](https://www.ugpti.org/resources/reports/downloads/mpc17-327.pdf)
  (report and paper share structure; sections below refer to the report).
- Redundancy = two complementary measures: **travel alternative diversity**
  and **network spare capacity** — report §2 (intro of section) and §1.2.
- **Effective route** = every link satisfies both:
  - *Efficient* (Dial 1971): moves strictly farther from origin,
    `l_r(head_a) > l_r(tail_a)` — report §2.1, Eq. (2-1).
  - *Not-too-long* (Leurent 1997): `(1+φ)(l_r(head_a) − l_r(tail_a)) ≥ l_a`;
    summing links bounds route length ≤ (1+φ)·shortest — report §2.1,
    Eqs. (2-2)–(2-3).
  - Elongation ratio **φ = 1.3–1.5 urban, 1.6 inter-urban**, citing
    Tagliacozzo & Pirzio (1973) and Leurent (1997) — report §2.1, text
    right after Eq. (2-2).
- Route counting: polynomial adjacency-matrix method of Meng et al. (2005)
  — report §2.1 and §1.2.
- Overlap treated only *indirectly* via link multiplicity `N_a^rs`
  (# effective routes through link a; concept from Russo & Vitetta 2003)
  — report §2.1, paragraph before Remark 1. **Gap we exploit: overlap is a
  diagnostic, not part of the index.**
- Spare capacity: reserve-capacity multiplier μ with logit/C-logit
  equilibrium — report §2.2. Needs OD demand (not demand-free).
- Case study: Winnipeg — report §1.2 (objectives) and §4.

### Rohrer, Jabbar & Sterbenz — Path diversification (EPD)
- [Telecommunication Systems 56:49–67 (2014)](https://link.springer.com/article/10.1007/s11235-013-9818-7);
  open PDF: [author copy](https://cdn.jprohrer.org/documents/publications/Rohrer-Jabbar-Sterbenz-2012.pdf).
- Path as element set `P = L ∪ N` (links AND nodes; node sharing matters —
  Baltimore tunnel fire example) — §3.1, Def. 1 and Fig. 1 discussion.
- **Path diversity** `D(P_b,P_a) = 1 − |P_b ∩ P_a|/|P_a|` — §3.1, Defs. 2–3,
  Eqs. (2)–(3).
- **Effective Path Diversity** `EPD = 1 − e^(−λ k_sd)`,
  `k_sd = Σ_i D_min(P_i)` (marginal diversity vs already-selected paths);
  λ experimentally set, λ=1 used — §3.2, Def. 3, Eqs. (4)–(5).
- **Total Graph Diversity** = mean EPD over node pairs; ring ≈ 0.6 at λ=1 —
  §3.5.
- Geographic diversity `D_g = α d²_min + βA` (min distance + enclosed area)
  — §3.3, Def. 4, Eq. (6). Barely developed; interesting for spatial nets.
- Validated only on internet topologies of **7–361 nodes** — §4, Table 1.
  **Gap we exploit: never applied to large planar street networks.**

## 2. Route-choice behavior (the empirical route envelope)

### Zhu & Levinson (2015) — Do people use the shortest path?
- [PLOS ONE 10(8):e0134322](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0134322).
- Data: GPS loggers (QSTARZ BT-Q1000p, point per 25 m), Minneapolis–St
  Paul; 143 subjects, **25,157 trips** (6,059 commutes), 8–13 weeks —
  Methods (data collection subsection).
- **~34% of all trips follow the shortest-time path** (exact match; ~40%
  with 10% overlap allowance); commutes only **13.5%** (21.7% at 10%
  allowance) — Results; also stated in Abstract.
- Time deviations are small: **50% of trips < 30 s** over shortest, **90%
  < 5 min**; non-commute 55% within +5% / 80% within +20% time; commute
  30% within +5% / 70% within +20% — Results (travel-time difference
  distributions).
- **Route portfolios: 249 distinct routes in 657 home–work trips from 95
  subjects (~2–3 routes per OD pair, used interchangeably)** — Results
  (three-week choice-set analysis).
- Deviation from shortest grows with trip length, overlap minimum ~20-mile
  trips — Results (deviation vs distance).
- Explanations: bounded rationality, reliability/road-type preferences,
  imperfect information, undetected side stops — Discussion.
- ⚠ caveat: one US metro, car commuters.

### FHWA multiday GPS study (Lexington, KY)
- [TMIP report, Chapter 4](https://www.fhwa.dot.gov/planning/tmip/publications/other_reports/multiday_gps/chapter04.cfm).
- <40% of commuters on shortest path; >50% of auto commutes exceed
  shortest time by 30%+; mean deviation 254–273 s; only ~25/124 trips
  with >90% overlap with shortest route — Chapter 4 (all).

## 3. Real traffic flows & demand structure

### Wang, Hidalgo et al. (2012) — Understanding road usage patterns
- [Scientific Reports 2:1001](https://www.nature.com/articles/srep01001);
  full text: [PMC3526957](https://pmc.ncbi.nlm.nih.gov/articles/PMC3526957/).
- Demand: phone billing records (360k users Bay Area = 6.56% pop; 680k
  Boston = 19.35%), trip = same user in two zones within 1 h; scaled by
  zone population share and vehicle-usage rates → transient OD (t-OD) —
  Methods (t-OD construction).
- **Routes are modeled, not observed**: Incremental Traffic Assignment,
  demand loaded 40/30/20/10%, travel times updated via BPR
  `t = t_f(1+0.15·VOC⁴)` — Methods (ITA subsection).
- Validation: predicted segment travel times vs GPS **probe vehicles**,
  R² > 0.9 — Methods (validation). Validates times, NOT route choices.
- Flow distribution P(V) ≈ **sum of two exponentials** (arterial/highway
  regimes; Bay Area means 373 vs 1,493 veh/h) — Results, Fig. 1.
- **98% of segments well below designed capacity** (VOC distribution) —
  Results, Fig. 1 (VOC panel).
- **Major driver sources**: each segment's 80%-of-flow catchment ≈ only
  ~20 zones (log-normal); each source touches ~1,000 segments (normal) —
  Results, "network of road usage", Fig. 2.
- Betweenness only weakly correlated with actual flow/K_road — Results
  (road-usage network section).

### Kirkley et al. (2018) — Betweenness invariance in street networks
- [Nature Communications 9:2501](https://www.nature.com/articles/s41467-018-04978-z).
- Betweenness distribution in 97 world cities is **bimodal & invariant**:
  a high-betweenness tree-like backbone carries the bulk of shortest
  paths; low-betweenness loop edges provide local alternatives — Results,
  Figs. 1–2. ⚠ section detail from abstract/figures, verify on full read.
- Explains why finite OD samples cover only a minority of edges
  (our 35–47% coverage observation).

## 4. Trip distances (locality)

- ⚠ Average urban trip ≈ 4.4 km; majority < 5 km; commute peak 2.5–7.5 km;
  trip lengths ~ truncated power law `P(l) ~ l^(−α)e^(−l/β)` — from
  [intra-urban mobility arXiv:1505.07372](https://arxiv.org/pdf/1505.07372),
  [urban scaling arXiv:1401.0207](https://arxiv.org/pdf/1401.0207),
  [tour behavior arXiv:2505.20590](https://arxiv.org/pdf/2505.20590).
  Search-level evidence; pin exact numbers before citing in the paper.

## 5. Reviews (framing / gap)

- ⚠ [Systematic analysis of road network vulnerability & resilience (2025)](https://www.sciencedirect.com/science/article/abs/pii/S2950024925000113),
  [SLR with Bibliometrix, 594 papers 2003–2023 (2025)](https://www.sciencedirect.com/science/article/pii/S2666188825007063),
  [Resilience to natural hazards SLR (2026)](https://www.sciencedirect.com/science/article/pii/S2666691X26000011)
  — metrics split topological/functional/hybrid; **limited link-level
  analysis flagged as gap**. From abstracts; verify claims on full read.
- Per-segment diversion-route index (closest per-edge relative):
  [Sustainability 14(4):2244 (2022)](https://doi.org/10.3390/su14042244)
  — vulnerability from two shortest disjoint diversion routes per segment.
  ⚠ Full text fetch blocked (403); read via journal site.

## 6. Our computational findings so far (2 km cores, 10 cities)

- Local directed endpoint index saturates in one-way systems: São Paulo
  72.5% of interior edges without directed alternative within +500 m; 32%
  of closures make u→v unreachable entirely — `src/redundancy.py` runs,
  2026-07-10.
- Trip-level index (800 uniform ODs, δ=1.3) restores discrimination:
  São Paulo mean V=0.605 (8.5% stranding edges), Tokyo 0.441 (3.1%) —
  `src/trip_index.py` runs, 2026-07-10.
- Coverage of uniform OD sampling: 35–47% of edges — consistent with
  backbone concentration (Kirkley et al.). Uncovered ≠ safe; motivates
  per-edge catchment sampling (v2).
- **Catchment-sampled index (v2, `src/catchment_index.py`, 2026-07-10)**:
  per-edge OD sampling from reverse/forward Dijkstra catchments; a trip is
  accepted when its through-edge route is within δ of its own shortest
  (tolerance envelope, not strict shortest-path membership); K=20 trips per
  edge, ≤2 per origin (Wang et al. ~20-zone catchment). **Coverage 35–47% →
  99.8–100%** (São Paulo 1,649/1,652 interior edges, Tokyo 2,769/2,770;
  no-demand edges 3 and 1).
- Unweighted trip aggregation compresses the index: marginal trips (shortest
  path survives closure) cap suffering at e^-λ ≈ 0.37 and dilute dependent
  trips — Tokyo 10–90 pct span was 0.29–0.35. Fix folded into v2: weight
  each trip by reliance `w = (δ − T/L0)/(δ − 1)` ∈ (0,1] (1 = through-edge
  route is the shortest; same linear-in-detour form as alternative quality).
- Weighted v2 results: São Paulo mean V=0.441, 10–90 pct 0.35–0.63,
  Spearman vs v1 rises 0.44→0.62 restricting to v1 edges with ≥10 trips
  (much v1–v2 disagreement is v1 sampling noise; v1 median trips/edge = 6).
  Tokyo mean V=0.327, span stays narrow (0.28–0.37) at any weighting —
  interpretable as real grid redundancy: mean reliance per accepted trip
  0.51 (Tokyo) vs 0.63 (São Paulo), i.e. in a grid most plausible users of
  an edge have near-equal alternate routes. Newly covered (off-backbone)
  edges score lower on average (SP 0.40) but include genuine stranding
  edges (~0.4%) — uncovered ≠ safe, confirmed.
- v2 runtime: single-pair alternative queries (`_shortest_alt`, 66% of
  runtime in a Tokyo profile) swapped from cutoff-Dijkstra to **budget-pruned
  A\*** with straight-line heuristic (`src/redundancy.py::_astar`,
  2026-07-10): admissible (edge length ≥ chord; equirectangular projection,
  0.999 safety factor) so results are exact — full-city CSVs byte-identical
  pre/post swap; 600-query harness 0 mismatches, 4.4× per query. The
  budget-pruning (`g + h > δ·L0` ⇒ prune) confines the search to the
  tolerance ellipse around (o, d), which is what pays in the frequent
  "no tolerable alternative" outcome. City wall-clock: São Paulo 194→118 s
  (72 ms/edge), Tokyo 878→425 s (154 ms/edge). Remaining bottleneck is now
  the multi-target catchment/origin Dijkstra sweeps (~⅓ of time pre-swap);
  candidate lever if more speed is needed: early-terminating sweeps or a
  compiled backend. Cutoff-aware bidirectional Dijkstra benchmarked as the
  alternative (identical results): statistical tie — São Paulo pipeline −12%,
  Tokyo +14% vs A\*, and A\* 1.5× ahead on a wide-budget random mix. One-way
  systems weaken the straight-line heuristic, grids strengthen it. A\* kept.
- Multi-target sweep lever taken next (2026-07-10): sweeps are keyed by
  *source node*, and interior nodes recur as edge endpoints / trip origins
  across many edges (Tokyo: ~28k origin sweeps over ~1.2k distinct sources,
  ~23× recomputation). `_SweepCache` in `src/catchment_index.py` caches each
  source's sweep once at maximal cutoff, storing results in Dijkstra settle
  order so smaller-radius node lists are exact prefixes — outputs (and the
  RNG stream) bit-identical to the naive implementation; full-city CSVs
  byte-equal. Rejected alternative: per-candidate A\* for L0 tests loses —
  Tokyo tolerance-tests ~18.6 candidates/origin at ~1.2 ms/query vs 6.6 ms
  for one amortized sweep. City wall-clock: São Paulo 118→58 s (35 ms/edge),
  Tokyo 425→195 s (70 ms/edge). Cumulative vs pre-optimization: SP 194→58 s
  (3.3×), Tokyo 878→195 s (4.5×). Remaining cost is now dominated by
  `trip_suffering` alternative extraction (~40–45%), already A\*-optimized —
  further speed means a compiled backend, not algorithmics.
- **Compiled backend + scale-up to 5 km cores (2026-07-10)**: to scale the
  study radius to the literature-grounded locality (~5 km trips), the v2
  pipeline was ported to a `scipy.sparse.csgraph` CSR core
  (`catchment_index.py::CSRCore`). Backend race on the real workloads
  (single-pair budget query with closure; cutoff sweep): scipy 0.26 / 0.31 ms
  vs igraph 0.40 / 1.05 ms (no native cutoff; mishandles inf weights) vs
  rustworkx 0.53 ms (Python cost callbacks) vs tuned-Python A\* 0.99 /
  nx 4.8 ms. Closures poke the collapsed CSR entry's weight; each entry
  keeps the multiset of open parallel lengths so per-key removal semantics
  match the MultiDiGraph exactly. Validation: 300 same-seed edges,
  **300/300 bit-exact V** vs the NetworkX implementation, 3.4× faster at
  2 km. Config bumped: download radius 6 km, buffer 1 km (scored interior
  = 5 km core), MAX_TRIP_M 3→5 km per the trip-length literature (§4).
  Caveat to state in the paper: alternatives budget δ·L0 ≤ 6.5 km can
  exceed the 1 km buffer for boundary-adjacent trips — cropping risk
  persists at the rim, as at 2 km.
- First 5 km-core results (2026-07-10): São Paulo 16,281/16,284 interior
  edges scored (100%, 3 no-demand), mean V=0.364, 16 min; Tokyo
  31,232/31,238 (100%, 6 no-demand), mean V=0.279, 64 min. Coverage holds
  at scale. Mean V drops vs the 2 km cores (SP 0.441→0.364, Tokyo
  0.327→0.279): longer trips (≤5 km) carry more absolute slack under the
  relative tolerance δ, and the periphery is less one-way-saturated than
  the historic cores — a scale effect to discuss alongside calibration.
  City contrast persists (SP mean and spread ~1.3–1.4× Tokyo's); mean
  reliance per trip now ≈0.58 in both (the 2 km SP/Tokyo gap 0.63/0.51
  was a core-area effect, not a city-wide one).

## Design decisions locked so far

| Parameter | Value | Anchor |
|---|---|---|
| Stretch bound δ | 1.2–1.3 (test 1.1–1.6) | Zhu & Levinson envelope; Xu & Chen φ |
| Route portfolio size | saturate after ~3 | Zhu & Levinson 249/657 routes |
| Quality weighting | decay with extra time | Zhu & Levinson 50% <30 s |
| Locality radius | ~5 km (2 km cores for now) | trip-length literature (⚠ pin) |
| Directionality | directed, always | our São Paulo saturation finding |
| Aggregation | EPD-style `1−e^(−λk)` | Rohrer & Sterbenz §3.2 |
| OD sampling | per-edge catchment, K≈20 trips, ≤2/origin | Wang et al. ~20-zone catchment; our coverage finding |
| Trip weighting | reliance `w=(δ−T/L0)/(δ−1)` | Zhu & Levinson tolerance envelope; our compression finding |
