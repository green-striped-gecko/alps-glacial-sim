"""
glacial_sim.py
==============

Spatially explicit coalescent simulation of Pleistocene glacial cycles in a
mountain landscape, used for the Australian Alps phylogeography paper.

The model is a 13 x 13 stepping-stone lattice of demes simulated backwards in
time with msprime. Every deme has its own effective size, which is rewritten
every 10,000 years according to a climate cycle and the deme's elevation. Demes
that fall outside the thermal tolerance of the modelled taxon become "ghost"
demes (Ne = 2, migration = 1e-5) so that lineages cannot pass through them, but
the population objects still exist and msprime stays happy.

Three landscape scenarios and two thermal ecologies are crossed:

  scenario        landscape
  --------        ---------
  "1peak"         one central massif, elevation falls off linearly with
                  distance from the grid centre
  "3peak"         three Gaussian massifs on an equilateral triangle
  "3peak_river"   as "3peak", plus an east-west river barrier that is active
                  during glacial phases

  species         thermal ecology
  -------         ---------------
  "warm"          lowland / warm-adapted: suitability increases with local
                  temperature; glacial refugia are the four grid corners
  "cold"          alpine / cold-adapted: suitability peaks at intermediate
                  local temperature; glacial refugia are the summits

Everything in this file is plain Python and can be run without R. The R layer
in ../R/ only wraps these functions through reticulate so that results come
back as genlight objects and data frames.

Time convention
---------------
All times are in generations before present, i.e. msprime's backwards-in-time
convention. t = 0 is the present day, which the model places at a warm peak
(interglacial).

Public entry points
-------------------
make_elevation(n, scenario)      elevation grid, values in [0, 1]
climate(t, ...)                  climate index at time t, 1 = warm, 0 = glacial
occupied_warm / occupied_cold    set of occupied deme indices at a given time
build_demography(...)            msprime.Demography for one scenario/species
sim_present(...)                 present-day genotype matrix
sim_timepoints(...)              genotype matrices at a series of time points
sim_pi(...)                      pi, He, Tajima's D, theta_W and Fst over time
demog_summary(...)               climate / Ne / occupancy bookkeeping, no msprime

Requires: numpy, msprime (>= 1.2), tskit.
"""

import math
import sys
import time

import numpy as np
import msprime

SCENARIOS = ("1peak", "3peak", "3peak_river")
SPECIES = ("warm", "cold")


# =============================================================================
# 1. LANDSCAPE
# =============================================================================

def _peak_summits(n):
    """Row/column indices of the three summits of the 3-peak landscapes.

    The summits sit on an equilateral triangle of radius n/3 centred on the
    grid, with one peak pointing north (90 degrees) and two to the south-west
    and south-east (210 and 330 degrees).
    """
    cx, cy = (n - 1) / 2.0, (n - 1) / 2.0
    r = n / 3.0
    peaks = []
    for angle_deg in (90, 210, 330):
        angle = math.radians(angle_deg)
        py = cy - r * math.sin(angle)
        px = cx + r * math.cos(angle)
        peaks.append((round(py), round(px)))
    return peaks


def make_elevation(n, scenario="1peak"):
    """Return an n x n elevation grid scaled to [0, 1].

    "1peak"        : elevation = 1 - d / d_max, where d is Euclidean distance
                     from the grid centre. A single cone, summit = 1 at the
                     centre, 0 at the corners.
    "3peak"        : maximum of three Gaussian bumps (sigma = n/5) placed at
                     _peak_summits(). Summits = 1, saddles between peaks stay
                     at intermediate elevation, so the three massifs remain
                     connected at moderate elevation.
    "3peak_river"  : identical elevation to "3peak"; the river is a barrier
                     applied to occupancy and migration, not to elevation.
    """
    _check_scenario(scenario)
    yy, xx = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")

    if scenario == "1peak":
        c = (n - 1) / 2.0
        d = np.sqrt((yy - c) ** 2 + (xx - c) ** 2)
        dmax = d.max() if d.max() > 0 else 1.0
        return 1.0 - d / dmax

    sigma = n / 5.0
    elev = np.zeros((n, n))
    for py, px in _peak_summits(n):
        d = np.sqrt((yy - py) ** 2 + (xx - px) ** 2)
        elev = np.maximum(elev, np.exp(-d ** 2 / (2 * sigma ** 2)))
    return elev


def _check_scenario(scenario):
    if scenario not in SCENARIOS:
        raise ValueError(
            "scenario must be one of %s, got %r" % (list(SCENARIOS), scenario))


# =============================================================================
# 2. CLIMATE
# =============================================================================

def climate(t, cycle_len=100_000, warm_len=20_000):
    """Asymmetric sawtooth climate index at time t (generations before present).

    Returns 1.0 at a warm peak (interglacial) and 0.0 at a glacial maximum.

    Read forwards in time the cycle is a slow 80 ky cooling into a glacial
    maximum followed by a fast 20 ky warming (`warm_len`) back to the
    interglacial. Because msprime runs backwards, t = 0 is the present
    interglacial and the first `warm_len` generations back from any warm peak
    are the fast warming limb traversed in reverse.
    """
    x = int(t) % cycle_len
    if x <= warm_len:
        return 1.0 - x / warm_len
    return (x - warm_len) / (cycle_len - warm_len)


# =============================================================================
# 3. THERMAL NICHE AND SUITABILITY
# =============================================================================
# Local temperature of a cell is  lt = climate - ELEV_WEIGHT * elevation,
# so high cells are colder than low cells at the same point in the cycle.
#
# Suitability is continuous in [0, 1] and maps linearly onto deme size:
#
#     Ne_cell = Ne_bad + (Ne_good - Ne_bad) * suitability
#
# so suitability 0 gives a refugium-sized deme and suitability 1 an
# optimum-sized deme, with no binary good/bad switch anywhere in the model.

ELEV_WEIGHT = 0.7      # how strongly elevation depresses local temperature
LT_MIN_WARM = 0.35     # warm-adapted: lower thermal limit for occupancy
ELEV_MAX_WARM = 0.70   # warm-adapted: upper elevation limit (below alpine zone)
LT_MIN_COLD = -0.55    # cold-adapted: lower thermal limit (not glaciated)
LT_MAX_COLD = 0.50     # cold-adapted: upper thermal limit (not too warm)
SUIT_HALF = 0.55       # cold-adapted: half-width of the suitability triangle

RIVER_THRESHOLD = 0.5  # river is active while climate <= this value


def suit_warm(elev_f, c):
    """Warm-adapted suitability: warmer is better, clipped to [0, 1]."""
    lt = c - ELEV_WEIGHT * elev_f
    return np.clip(lt, 0.0, 1.0)


def suit_cold(elev_f, c):
    """Cold-adapted suitability: a tent function peaking at lt = 0.

    Suitability falls to 0 at lt <= -SUIT_HALF (glaciated) and at
    lt >= +SUIT_HALF (too warm). At an interglacial summit lt is about 0.30,
    giving suitability about 0.45, i.e. refugium-level deme sizes; during a
    mild glacial the summit sits near lt = 0 and reaches full size.
    """
    lt = c - ELEV_WEIGHT * elev_f
    return np.clip(1.0 - np.abs(lt) / SUIT_HALF, 0.0, 1.0)


# =============================================================================
# 4. RIVER BARRIER  (scenario "3peak_river" only)
# =============================================================================

def river_row(n):
    """Grid row occupied by the river, 0-based: (n // 2) - 1.

    For n = 13 this is row 5, i.e. one row north of the grid centre, which
    separates the northern massif from the two southern massifs.
    """
    return (n // 2) - 1


def _apply_river(occ, n, c, scenario):
    """Remove the river row from the occupied set while the river is active.

    The river is glacially fed, so it only exists while climate <=
    RIVER_THRESHOLD. Cells on the river row are then forced to ghost status.
    This is applied *after* refugia are force-included, so a refugium that
    happens to sit on the river row is still removed.
    """
    if scenario == "3peak_river" and c <= RIVER_THRESHOLD:
        rr = river_row(n)
        for cj in range(n):
            occ.discard(rr * n + cj)
    return occ


# =============================================================================
# 5. OCCUPANCY
# =============================================================================

def occupied_warm(elev, c, scenario="1peak"):
    """Deme indices occupied by the warm-adapted taxon at climate c.

    A cell is occupied if it is warm enough (lt >= LT_MIN_WARM) and low enough
    (elevation <= ELEV_MAX_WARM). During glacial phases (c <= 0.5) the four
    corner refugia and their four-neighbours are force-included so that
    lineages survive the glacial maximum in the lowland corners.
    """
    _check_scenario(scenario)
    n = elev.shape[0]
    lt = c - ELEV_WEIGHT * elev
    elev_f = elev.flatten()
    mask = (lt.flatten() >= LT_MIN_WARM) & (elev_f <= ELEV_MAX_WARM)
    occ = set(int(k) for k in np.where(mask)[0]) if mask.any() else set()

    if c <= 0.5:
        for ci, cj in [(0, 0), (0, n - 1), (n - 1, 0), (n - 1, n - 1)]:
            occ.add(ci * n + cj)
            for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                ni, nj = ci + di, cj + dj
                if 0 <= ni < n and 0 <= nj < n:
                    occ.add(ni * n + nj)

    return _apply_river(occ, n, c, scenario)


def _cold_refugia_cells(n, scenario):
    """Summit refugia for the cold-adapted taxon.

    For each summit (one for "1peak", three otherwise) this is the summit cell,
    its four-neighbour ring, and the four cells at orthogonal distance two.
    Duplicates are removed because neighbouring massifs can share cells.
    """
    _check_scenario(scenario)
    offsets = [(0, 0),
               (-1, 0), (1, 0), (0, -1), (0, 1),
               (-2, 0), (2, 0), (0, -2), (0, 2)]
    summits = [(n // 2, n // 2)] if scenario == "1peak" else _peak_summits(n)

    cells = []
    for ci, cj in summits:
        for di, dj in offsets:
            r, c_ = ci + di, cj + dj
            if 0 <= r < n and 0 <= c_ < n:
                cells.append(r * n + c_)
    return list(set(cells))


def occupied_cold(elev, c, t_abs=None, cycle_len=100_000, scenario="1peak"):
    """Deme indices occupied by the cold-adapted taxon at climate c.

    A cell is occupied if its local temperature lies inside the thermal window
    [LT_MIN_COLD, LT_MAX_COLD]: not glaciated, not too warm.

    Near the glacial maximum the summit refugia are force-included even if the
    thermal mask would drop them, so that the alpine lineages persist. The
    window is t_abs mod cycle_len in [10 ky, 40 ky] combined with c <= 0.25
    (equivalently 60-90 ky into the cycle read forwards), plus an
    unconditional force at the deepest glacial (c <= 0.1). Refugia are only
    added if they are not glaciated (lt >= LT_MIN_COLD).
    """
    _check_scenario(scenario)
    n = elev.shape[0]
    lt = c - ELEV_WEIGHT * elev
    lt_f = lt.flatten()
    mask = (lt_f >= LT_MIN_COLD) & (lt_f <= LT_MAX_COLD)
    occ = set(int(k) for k in np.where(mask)[0]) if mask.any() else set()

    in_lgm_window = False
    if t_abs is not None:
        t_mod = int(t_abs) % cycle_len
        in_lgm_window = (10_000 <= t_mod <= 40_000) and (c <= 0.25)

    if in_lgm_window or c <= 0.1:
        for k in _cold_refugia_cells(n, scenario):
            if lt_f[k] >= LT_MIN_COLD:
                occ.add(k)

    return _apply_river(occ, n, c, scenario)


def occupied(elev, c, species, scenario, t_abs=None, cycle_len=100_000):
    """Dispatch to occupied_warm or occupied_cold."""
    if species == "warm":
        return occupied_warm(elev, c, scenario)
    return occupied_cold(elev, c, t_abs=t_abs, cycle_len=cycle_len,
                         scenario=scenario)


# =============================================================================
# 6. GRID HELPERS
# =============================================================================

def pop_name(k):
    """msprime population name for deme index k."""
    return f"p{k}"


def cell_idx(i, j, n):
    """Row-major deme index of grid cell (row i, column j)."""
    return i * n + j


def neighbors4(k, n):
    """Yield the four-neighbour deme indices of deme k on an n x n grid."""
    i, j = divmod(k, n)
    for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        ni, nj = i + di, j + dj
        if 0 <= ni < n and 0 <= nj < n:
            yield cell_idx(ni, nj, n)


def time_grid(total, dt=10_000):
    """Times at which the demography is rewritten: every dt up to total."""
    return list(range(0, total + 1, dt))


# =============================================================================
# 7. DEME SIZES
# =============================================================================

NE_GHOST = 2       # size of an unoccupied deme
M_GHOST = 1e-5     # migration into or out of an unoccupied deme


def cell_ne_grid(species, elev_f, c, occ, Ne_good, Ne_bad, Ne_ghost=NE_GHOST):
    """Vector of deme sizes for all n*n cells at climate c.

    Occupied cells get Ne_bad + (Ne_good - Ne_bad) * suitability, floored at
    Ne_ghost. Unoccupied cells stay at Ne_ghost.
    """
    s = suit_warm(elev_f, c) if species == "warm" else suit_cold(elev_f, c)
    ne = np.full(len(elev_f), float(Ne_ghost))
    for k in occ:
        ne[k] = max(Ne_bad + (Ne_good - Ne_bad) * s[k], float(Ne_ghost))
    return ne


# =============================================================================
# 8. DEMOGRAPHY
# =============================================================================

def build_demography(n, species, total_time,
                     Ne_good, Ne_bad,
                     m_connected, m_isolated,
                     cycle_len, warm_len,
                     scenario="1peak", t_offset=0):
    """Build the full msprime.Demography for one scenario and species.

    One population per grid cell. Every 10,000 generations (see time_grid) the
    size of each deme and the migration rate on each lattice edge are rewritten
    to match the climate at that time.

    Migration rates on an edge between two occupied demes:
      m_connected   default
      m_isolated    warm-adapted taxon during a glacial (c <= 0.35), when the
                    surviving demes are the four isolated corner refugia
      m_isolated    "3peak_river" only: any edge touching the river row while
                    the river is active
    Edges touching an unoccupied deme always get M_GHOST.

    t_offset shifts the climate clock, so that a simulation whose samples are
    taken at time T still sees the correct climate history: pass t_offset = T.

    Returns (demography, elevation_grid).
    """
    _check_scenario(scenario)
    elev = make_elevation(n, scenario)
    elev_f = elev.flatten()
    npops = n * n

    dem = msprime.Demography()
    for k in range(npops):
        dem.add_population(name=pop_name(k), initial_size=NE_GHOST)
    for k in range(npops):
        for nb in neighbors4(k, n):
            dem.set_migration_rate(pop_name(k), pop_name(nb), 0.0)

    rr = river_row(n)

    for t in time_grid(total_time):
        t_abs = t + t_offset
        c = climate(t_abs, cycle_len, warm_len)
        occ = occupied(elev, c, species, scenario,
                       t_abs=t_abs, cycle_len=cycle_len)
        ne = cell_ne_grid(species, elev_f, c, occ, Ne_good, Ne_bad)

        warm_isolated = (species == "warm" and c <= 0.35)
        river_active = (scenario == "3peak_river" and c <= RIVER_THRESHOLD)

        for k in range(npops):
            dem.add_population_parameters_change(
                time=t, population=pop_name(k),
                initial_size=int(max(ne[k], 1)))

        for k in range(npops):
            for nb in neighbors4(k, n):
                if (k in occ) and (nb in occ):
                    touches_river = river_active and (
                        k // n == rr or nb // n == rr)
                    if touches_river:
                        rate = m_isolated
                    else:
                        rate = m_isolated if warm_isolated else m_connected
                else:
                    rate = M_GHOST
                dem.add_migration_rate_change(
                    time=t, source=pop_name(k), dest=pop_name(nb), rate=rate)

    return dem, elev


def _burn_in(n, Ne_good, cycle_len):
    """Extra time appended after the modelled cycles so lineages can coalesce.

    Rounded up to a whole number of climate cycles from a rough expectation of
    the time to the most recent common ancestor in a structured population,
    4 * (2 * Ne_good * n * n).
    """
    t_mrca_est = 2 * Ne_good * n * n
    return int(np.ceil(4 * t_mrca_est / cycle_len)) * cycle_len


def _sample_sets(occ, samples_per_cell):
    """One haploid SampleSet of 2 * samples_per_cell nodes per occupied deme.

    Nodes are drawn as haploids and paired into diploid genotypes afterwards,
    which keeps the mapping from sample nodes to demes explicit.
    """
    return [msprime.SampleSet(num_samples=2 * samples_per_cell,
                              population=pop_name(k), time=0, ploidy=1)
            for k in sorted(occ)]


# =============================================================================
# 9A. PRESENT-DAY GENOTYPES
# =============================================================================

def sim_present(n, species, total_time,
                samples_per_cell, seq_len, mu, seed,
                Ne_good, Ne_bad, m_connected, m_isolated,
                cycle_len, warm_len, scenario="1peak"):
    """Simulate present-day (t = 0) genotypes for every occupied deme.

    Returns a dict with
      G        integer genotype matrix, individuals x sites, values 0/1/2
      pos      site positions along the simulated sequence
      ind_pop  deme name for each individual
      n_ind, n_sites
    """
    dem, elev = build_demography(n, species, total_time,
                                 Ne_good, Ne_bad,
                                 m_connected, m_isolated,
                                 cycle_len, warm_len,
                                 scenario=scenario, t_offset=0)

    c0 = climate(0, cycle_len, warm_len)
    occ = occupied(elev, c0, species, scenario, t_abs=0, cycle_len=cycle_len)

    ind_pop = []
    for k in sorted(occ):
        ind_pop.extend([pop_name(k)] * samples_per_cell)

    ts = msprime.sim_ancestry(samples=_sample_sets(occ, samples_per_cell),
                              demography=dem,
                              sequence_length=seq_len,
                              recombination_rate=1e-8,
                              random_seed=seed)
    mts = msprime.sim_mutations(ts, rate=mu,
                                model=msprime.BinaryMutationModel(),
                                random_seed=seed + 1)

    G = mts.genotype_matrix().astype(np.int8)
    G_ind = (G[:, 0::2] + G[:, 1::2]).T
    pos = np.array([v.position for v in mts.variants()])

    return {"G": G_ind, "pos": pos, "ind_pop": ind_pop,
            "n_sites": G_ind.shape[1], "n_ind": G_ind.shape[0]}


# =============================================================================
# 9B. GENOTYPES AT A SERIES OF TIME POINTS
# =============================================================================

def sim_timepoints(n, species, total_time,
                   sample_times, samples_per_cell,
                   seq_len, mu, seed,
                   Ne_good, Ne_bad, m_connected, m_isolated,
                   cycle_len, warm_len, scenario="1peak"):
    """Simulate genotypes for a series of sampling times.

    Each sampling time is an independent simulation whose climate clock is
    offset by that time, with the demography extended by a burn-in so that
    lineages have time to coalesce. Returns a list of per-deme blocks, one
    dict per (time, deme).

    This is switched off in the published runs (do_gl_timeseries = FALSE);
    sim_pi below computes the summary statistics directly from the tree
    sequences instead, which is far cheaper.
    """
    _check_scenario(scenario)
    elev = make_elevation(n, scenario)
    results = []
    burn_in = _burn_in(n, Ne_good, cycle_len)

    n_tp = len(sample_times)
    t0_all = time.time()
    for i, t_sample in enumerate(sample_times):
        t0 = time.time()
        t_extra = total_time - int(t_sample) + burn_in

        dem, _ = build_demography(n, species, t_extra,
                                  Ne_good, Ne_bad, m_connected, m_isolated,
                                  cycle_len, warm_len,
                                  scenario=scenario, t_offset=int(t_sample))

        c = climate(int(t_sample), cycle_len, warm_len)
        occ = occupied(elev, c, species, scenario,
                       t_abs=int(t_sample), cycle_len=cycle_len)

        ts = msprime.sim_ancestry(samples=_sample_sets(occ, samples_per_cell),
                                  demography=dem,
                                  sequence_length=seq_len,
                                  recombination_rate=1e-8,
                                  random_seed=seed + i)
        mts = msprime.sim_mutations(ts, rate=mu,
                                    model=msprime.BinaryMutationModel(),
                                    random_seed=seed + i + 1)

        G = mts.genotype_matrix().astype(np.int8)
        pos = np.array([v.position for v in mts.variants()])

        col = 0
        for k in sorted(occ):
            ncols = 2 * samples_per_cell
            haps = G[:, col:col + ncols]
            G_blk = (haps[:, 0::2] + haps[:, 1::2]).T
            results.append({"t": t_sample, "pop": pop_name(k),
                            "G": G_blk, "pos": pos,
                            "n_ind": G_blk.shape[0],
                            "n_sites": G_blk.shape[1]})
            col += ncols

        _progress("sim_tp", species, i, n_tp, t_sample, len(occ), t0, t0_all)

    return results


# =============================================================================
# 9C. SUMMARY STATISTICS THROUGH TIME
# =============================================================================

def sim_pi(n, species, total_time,
           sample_times, samples_per_cell,
           seq_len, mu, seed,
           Ne_good, Ne_bad, m_connected, m_isolated,
           cycle_len, warm_len, scenario="1peak"):
    """Nucleotide diversity and friends at each sampling time.

    For every time point an independent tree sequence is simulated exactly as
    in sim_timepoints, then summarised without ever materialising a genlight:

      pi        mean per-deme nucleotide diversity (tskit branch-free,
                site-based diversity, per base pair)
      tajima_d  mean per-deme Tajima's D
      theta_w   mean per-deme Watterson's theta, S / a1, in units of segregating
                sites (not divided by sequence length)
      he_mean   mean per-deme expected heterozygosity over segregating sites
      fst_mean  mean Hudson Fst over *four-neighbour deme pairs only*
                (adjacent pairs, not all pairs; the all-pairs present-day Fst
                used for H2 and H3 comes from R/05_fst_extraction.R)

    Demes with fewer than four sample nodes are skipped. A time point with no
    usable deme returns NaN for every statistic.
    """
    _check_scenario(scenario)
    elev = make_elevation(n, scenario)
    pis, taj_ds, theta_ws, fst_means, he_means = [], [], [], [], []
    burn_in = _burn_in(n, Ne_good, cycle_len)

    n_tp = len(sample_times)
    t0_all = time.time()
    for i, t_sample in enumerate(sample_times):
        t0 = time.time()
        t_extra = total_time - int(t_sample) + burn_in

        dem, _ = build_demography(n, species, t_extra,
                                  Ne_good, Ne_bad, m_connected, m_isolated,
                                  cycle_len, warm_len,
                                  scenario=scenario, t_offset=int(t_sample))

        c = climate(int(t_sample), cycle_len, warm_len)
        occ = occupied(elev, c, species, scenario,
                       t_abs=int(t_sample), cycle_len=cycle_len)

        ts = msprime.sim_ancestry(samples=_sample_sets(occ, samples_per_cell),
                                  demography=dem,
                                  sequence_length=seq_len,
                                  recombination_rate=1e-8,
                                  random_seed=seed + i)
        mts = msprime.sim_mutations(ts, rate=mu,
                                    model=msprime.BinaryMutationModel(),
                                    random_seed=seed + i + 1)

        # map sample nodes back to demes, in the order the sample sets were built
        all_nodes = list(mts.samples())
        pop_nodes = {}
        col = 0
        for k in sorted(occ):
            n_haps = samples_per_cell * 2
            pop_nodes[k] = all_nodes[col:col + n_haps]
            col += n_haps

        cell_node_lists = [v for v in pop_nodes.values() if len(v) >= 4]

        if not cell_node_lists:
            for lst in (pis, taj_ds, theta_ws, fst_means, he_means):
                lst.append(float("nan"))
            continue

        try:
            G_full = mts.genotype_matrix()
            node_to_col = {nd: ii for ii, nd in enumerate(all_nodes)}
        except Exception:
            G_full = None
            node_to_col = {}

        cell_pis, cell_taj, cell_thw, cell_he = [], [], [], []
        for nodes in cell_node_lists:
            cell_pis.append(float(mts.diversity(sample_sets=[nodes])[0]))
            cell_taj.append(float(mts.Tajimas_D(sample_sets=[nodes])[0]))
            s = float(mts.segregating_sites(sample_sets=[nodes])[0])
            a1 = sum(1.0 / j for j in range(1, len(nodes)))
            cell_thw.append(s / a1 if a1 > 0 else float("nan"))

            n_h = len(nodes)
            idx = [node_to_col[nd] for nd in nodes]
            if G_full is not None and G_full.shape[1] > 0:
                ac = G_full[:, idx].sum(axis=1).astype(float)
                p = ac / n_h
                seg = (p > 0) & (p < 1)
                he_seg = 2 * p[seg] * (1 - p[seg])
                cell_he.append(float(he_seg.mean()) if seg.any() else float("nan"))
            else:
                cell_he.append(float("nan"))

        pis.append(float(np.nanmean(cell_pis)))
        taj_ds.append(float(np.nanmean(cell_taj)))
        theta_ws.append(float(np.nanmean(cell_thw)))
        he_means.append(float(np.nanmean(cell_he)))

        _progress("sim_pi", species, i, n_tp, t_sample, len(occ), t0, t0_all)

        # Fst between four-neighbour deme pairs
        if len(cell_node_lists) >= 2:
            fst_vals = []
            for k in sorted(pop_nodes.keys()):
                for nb in neighbors4(k, n):
                    if nb > k and nb in pop_nodes:
                        na_, nb_ = pop_nodes[k], pop_nodes[nb]
                        if len(na_) >= 4 and len(nb_) >= 4:
                            try:
                                f = mts.Fst(sample_sets=[na_, nb_])
                                v = float(np.asarray(f).flat[0])
                                if not np.isnan(v):
                                    fst_vals.append(max(0.0, v))
                            except Exception:
                                pass
            fst_means.append(float(np.nanmean(fst_vals)) if fst_vals
                             else float("nan"))
        else:
            fst_means.append(float("nan"))

    return {"times": np.array(sample_times, dtype=int),
            "pi": np.array(pis, dtype=float),
            "tajima_d": np.array(taj_ds, dtype=float),
            "theta_w": np.array(theta_ws, dtype=float),
            "fst_mean": np.array(fst_means, dtype=float),
            "he_mean": np.array(he_means, dtype=float)}


def _progress(tag, species, i, n_tp, t_sample, n_occ, t0, t0_all):
    """One-line progress report on stderr with a running estimate of time left."""
    elapsed = time.time() - t0
    total_e = time.time() - t0_all
    eta = total_e / (i + 1) * (n_tp - i - 1)
    print(f"  {tag} [{species}] {i+1:>2}/{n_tp}  t={t_sample//1000}k  "
          f"n_occ={n_occ:>2}  {elapsed:.1f}s  ETA {eta/60:.1f}min",
          file=sys.stderr, flush=True)


# =============================================================================
# 10. DEMOGRAPHY BOOKKEEPING  (no coalescent simulation)
# =============================================================================

def demog_summary(n, species, total_time, times,
                  Ne_good, Ne_bad, cycle_len, warm_len, scenario="1peak"):
    """Climate, occupancy and deme sizes over time, without running msprime.

    Cheap enough to run on a fine time grid. Used for the landscape/occupancy
    figure and for sanity-checking the demography that the coalescent sees.

    Returns climate, summed Ne over occupied demes, number of occupied demes,
    the full time x deme size matrix, and the flattened elevation grid.
    """
    _check_scenario(scenario)
    elev = make_elevation(n, scenario)
    elev_f = elev.flatten()
    times = list(times)
    clim_arr = np.zeros(len(times))
    Ne_arr = np.zeros(len(times))
    n_occ_arr = np.zeros(len(times), dtype=int)
    Ne_cell = np.zeros((len(times), n * n))

    for it, t in enumerate(times):
        c = climate(t, cycle_len, warm_len)
        clim_arr[it] = c
        occ = occupied(elev, c, species, scenario, t_abs=t, cycle_len=cycle_len)
        ne = cell_ne_grid(species, elev_f, c, occ, Ne_good, Ne_bad)

        n_occ_arr[it] = len(occ)
        Ne_arr[it] = ne[list(occ)].sum() if occ else 0.0
        Ne_cell[it, :] = ne

    return {"times": np.array(times, dtype=int),
            "climate": clim_arr,
            "Ne_global": Ne_arr,
            "n_occ": n_occ_arr,
            "Ne_cell": Ne_cell,
            "elev": elev_f}
