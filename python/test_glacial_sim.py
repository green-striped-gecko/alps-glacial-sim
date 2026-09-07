"""
Sanity checks for glacial_sim.py.

Not a full test suite: these are the checks that catch the mistakes that
actually matter here, namely a landscape that does not look like the intended
landscape, a climate cycle that is out of phase, refugia that vanish at the
glacial maximum, and a river that fails to separate the northern massif from
the southern ones.

Run from the repository root:
    python python/test_glacial_sim.py
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import glacial_sim as gs

N = 13
CYCLE, WARM = 100_000, 20_000
failures = []


def check(label, cond):
    print(("  PASS  " if cond else "  FAIL  ") + label)
    if not cond:
        failures.append(label)


print("landscape")
for sc in gs.SCENARIOS:
    e = gs.make_elevation(N, sc)
    check(f"{sc}: shape and range",
          e.shape == (N, N) and 0.0 <= e.min() and e.max() <= 1.0 + 1e-12)

e1 = gs.make_elevation(N, "1peak")
check("1peak: summit at the grid centre",
      np.unravel_index(np.argmax(e1), e1.shape) == (N // 2, N // 2))

e3 = gs.make_elevation(N, "3peak")
check("3peak: three cells at elevation 1", np.isclose(e3, 1.0).sum() == 3)
check("3peak_river: same elevation as 3peak",
      np.allclose(e3, gs.make_elevation(N, "3peak_river")))

print("\nclimate")
check("interglacial at t = 0", gs.climate(0, CYCLE, WARM) == 1.0)
check("glacial maximum one warm_len back",
      gs.climate(WARM, CYCLE, WARM) == 0.0)
check("back to interglacial after one cycle",
      gs.climate(CYCLE, CYCLE, WARM) == 1.0)
c = [gs.climate(t, CYCLE, WARM) for t in range(0, 3 * CYCLE, 500)]
check("stays within [0, 1]", min(c) >= 0.0 and max(c) <= 1.0)

print("\noccupancy")
for sc in gs.SCENARIOS:
    elev = gs.make_elevation(N, sc)
    for sp in gs.SPECIES:
        occ_sizes = []
        for t in range(0, 3 * CYCLE, 2500):
            cc = gs.climate(t, CYCLE, WARM)
            occ = gs.occupied(elev, cc, sp, sc, t_abs=t, cycle_len=CYCLE)
            occ_sizes.append(len(occ))
        check(f"{sc}/{sp}: never empty", min(occ_sizes) > 0)
        check(f"{sc}/{sp}: contracts and expands",
              min(occ_sizes) < max(occ_sizes))

print("\nriver barrier")
elev = gs.make_elevation(N, "3peak_river")
rr = gs.river_row(N)
c_gl = 0.0
occ_gl = gs.occupied(elev, c_gl, "cold", "3peak_river",
                     t_abs=WARM, cycle_len=CYCLE)
check("river row empty at the glacial maximum",
      not any(k // N == rr for k in occ_gl))
# the warm-adapted taxon is the one that occupies the lowland river row at an
# interglacial; the cold-adapted taxon is confined to high ground and is never
# on that row, river or no river
occ_ig = gs.occupied(elev, 1.0, "warm", "3peak_river", t_abs=0, cycle_len=CYCLE)
check("river row occupied at the interglacial",
      any(k // N == rr for k in occ_ig))
occ_no = gs.occupied(elev, c_gl, "cold", "3peak", t_abs=WARM, cycle_len=CYCLE)
check("3peak keeps the row the river would remove", occ_no != occ_gl)

print("\ndeme sizes")
elev_f = gs.make_elevation(N, "1peak").flatten()
occ = gs.occupied(gs.make_elevation(N, "1peak"), 1.0, "warm", "1peak")
ne = gs.cell_ne_grid("warm", elev_f, 1.0, occ, 500, 50)
check("occupied demes between Ne_bad and Ne_good",
      all(50 - 1e-9 <= ne[k] <= 500 + 1e-9 for k in occ))
check("unoccupied demes are ghosts",
      all(ne[k] == gs.NE_GHOST for k in range(N * N) if k not in occ))

print("\nshort end-to-end simulation (5 x 5 grid, one cycle)")
res = gs.sim_pi(n=5, species="cold", total_time=CYCLE,
                sample_times=[0, 30_000, 60_000], samples_per_cell=4,
                seq_len=2e5, mu=1e-8, seed=1,
                Ne_good=100, Ne_bad=20, m_connected=0.005, m_isolated=0.0005,
                cycle_len=CYCLE, warm_len=WARM, scenario="3peak_river")
check("three time points returned", len(res["pi"]) == 3)
check("diversity is finite and positive", np.all(np.isfinite(res["pi"])) and
      np.all(res["pi"] > 0))
# fst_mean is NaN wherever a time point has fewer than two adjacent demes with
# enough samples, which happens on a grid this small; only finite values are
# constrained
fin = res["fst_mean"][np.isfinite(res["fst_mean"])]
check("finite Fst values lie in [0, 1]",
      len(fin) > 0 and np.all((fin >= 0) & (fin <= 1)))

pres = gs.sim_present(n=5, species="warm", total_time=CYCLE,
                      samples_per_cell=3, seq_len=2e5, mu=1e-8, seed=2,
                      Ne_good=100, Ne_bad=20, m_connected=0.005,
                      m_isolated=0.0005, cycle_len=CYCLE, warm_len=WARM,
                      scenario="1peak")
check("genotypes are 0/1/2", set(np.unique(pres["G"])) <= {0, 1, 2})
check("one population label per individual",
      len(pres["ind_pop"]) == pres["n_ind"])

check("same seed gives the same result",
      np.array_equal(pres["G"],
                     gs.sim_present(n=5, species="warm", total_time=CYCLE,
                                    samples_per_cell=3, seq_len=2e5, mu=1e-8,
                                    seed=2, Ne_good=100, Ne_bad=20,
                                    m_connected=0.005, m_isolated=0.0005,
                                    cycle_len=CYCLE, warm_len=WARM,
                                    scenario="1peak")["G"]))

print("\n%d failure(s)" % len(failures))
sys.exit(1 if failures else 0)
