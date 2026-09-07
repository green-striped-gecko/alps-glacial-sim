# Glacial-cycle simulations for the Australian Alps

Spatially explicit coalescent simulations used to ask how repeated Pleistocene
glacial cycles shape genetic diversity and population structure in a mountain
landscape, and how the answer depends on a taxon's thermal ecology.

Everything needed to reproduce the simulations and the simulation figures is
here. See **[`doc/SOFTWARE_REQUIREMENTS.md`](doc/SOFTWARE_REQUIREMENTS.md)** for what to
install and how (also available as
`doc/software_requirements.pdf`).

[![DOI](https://zenodo.org/badge/1359648290.svg)](https://doi.org/10.5281/zenodo.22581727)
---

## The model in one page

A **13 x 13 lattice of demes** is simulated backwards in time with
[msprime](https://tskit.dev/msprime/). Each cell of the grid has an elevation,
and a climate index oscillates through **five 100,000-generation glacial
cycles** (500 ky in total). The cycle is asymmetric: read forwards, an 80 ky
slow slide into a glacial maximum followed by a fast 20 ky warming back to an
interglacial. Time zero, the present day, sits at an interglacial peak.

Each cell's **local temperature** is the climate index minus 0.7 times its
elevation, so high cells are colder. Local temperature drives a continuous
**suitability** in [0, 1], which maps linearly onto the deme's effective size:

```
Ne_cell = Ne_bad + (Ne_good - Ne_bad) * suitability
```

with `Ne_good = 500` and `Ne_bad = 50`. There is no binary good/bad switch
anywhere in the model. Cells outside the taxon's thermal tolerance become
**ghost demes** (Ne = 2, migration 1e-5): the msprime population still exists,
so the demography stays valid, but lineages effectively cannot pass through.
Deme sizes and all migration rates are rewritten every 10,000 generations.

### Two thermal ecologies

| | Suitability | Glacial refugia |
|---|---|---|
| **warm-adapted** | rises with local temperature; lowland is best | the four lowland grid corners, which are mutually isolated (`m_isolated`) |
| **cold-adapted** | peaks at intermediate local temperature; too warm and too cold are both bad | the summits, which stay connected through saddle cells (`m_connected`) |

### Three landscapes

| Scenario | Landscape |
|---|---|
| `1peak` | one central massif; elevation falls linearly with distance from the centre |
| `3peak` | three Gaussian massifs on an equilateral triangle, connected by saddles |
| `3peak_river` | as `3peak`, plus an east-west river on grid row `(n %/% 2) - 1` that exists only while climate <= 0.5, i.e. during glacials. River-row cells are forced to ghost status and edges touching them drop to `m_isolated` |

The river separates the northern massif from the two southern ones. Note that
during the glacial phases in which the river exists, the warm-adapted taxon is
already confined to the lowland corners, so the barrier never touches it: the
warm-adapted arm of `3peak_river` is identical to `3peak` by construction, and
Figure 2 omits it by default.

### The design

3 landscapes x 30 replicates = **90 runs**, each simulating both thermal
ecologies, so 180 simulated data sets. Model parameters are identical
everywhere; only the landscape and the seed change. The same 30 seeds are
reused across scenarios, so replicate *r* of `1peak` and replicate *r* of
`3peak` are a matched pair. The exact seeds are in `R/02_design.R` and in
`design/published_design_matrix.csv`.

---

## Layout

```
python/
  glacial_sim.py            the simulation engine; pure Python, runs standalone
  test_glacial_sim.py       quick sanity checks on the engine

R/
  00_setup.R                repository paths, conda environment, version logging
  01_simulate.R             reticulate wrappers; run_simulation(), run_both()
  02_design.R               the design matrix and the published seeds
  03_run_batch.R            serial runner, load_batch(), extract_stat_matrix()
  04_run_parallel.R         parallel runner (PSOCK cluster)
  05_extract_fst.R          present-day all-pairs Hudson Fst and per-deme pi
  06_fig_h1.R               Figure 2  (H1, thermal ecology)
  07_fig_h2_h3.R            Figures 5 and 6  (H2 terrain, H3 river)
  08_fig1_occupancy.R       Figure 1  (occupancy through a glacial cycle)
  09_fig_supplementary.R    earlier scenario-comparison figures

scripts/
  00_check_install.R        verify the installation with a tiny test run
  01_run_simulations.R      run the design
  02_extract_summaries.R    build the pairwise / per-deme tables
  03_make_figures.R         build the figures

design/published_design_matrix.csv
doc/software_requirements.pdf
environment.yml
```

---

## Quick start

```bash
conda env create -f environment.yml
conda activate msprime-env
Rscript scripts/00_check_install.R        # tiny end-to-end test, under a minute
```

One version pairing matters: reticulate must be 1.38 or newer to hand numpy 2
arrays back to R, and msprime 1.4 onwards requires numpy 2. `init_python()`
checks this and stops with an explanation if the combination is wrong. See
`doc/SOFTWARE_REQUIREMENTS.md`.

Then, from the repository root:

```bash
Rscript scripts/01_run_simulations.R      # -> results/batch/out_*.rds
Rscript scripts/02_extract_summaries.R    # -> results/pw_dv.rds
Rscript scripts/03_make_figures.R         # -> figures/*.png, *.pdf
```

Figure 1 is the exception: it is drawn entirely from climate, occupancy and
per-cell Ne, so `build_occupancy_data()` regenerates its input in seconds and
it can be rebuilt without running any simulations at all.

Edit the settings block at the top of `scripts/01_run_simulations.R` before the
first real run: conda environment, output directory, number of replicates, number
of cores. Start with `N_REPS <- 1`.

To run a single simulation interactively:

```r
source("R/00_setup.R"); source("R/01_simulate.R")
init_python()

res <- run_simulation(species = "cold", scenario = "3peak_river",
                      n = 13, cycles = 5, seed = 1111)

head(res$pi)          # diversity, Tajima's D, theta_W, Fst through time
res$gl_present        # genlight: 20 diploids from every occupied deme at t = 0
```

The Python engine also runs on its own, without R:

```python
import sys; sys.path.insert(0, "python")
import glacial_sim as gs

out = gs.sim_pi(n=13, species="cold", total_time=500_000,
                sample_times=list(range(0, 500_000, 5000)),
                samples_per_cell=20, seq_len=2e6, mu=1e-8, seed=1311,
                Ne_good=500, Ne_bad=50, m_connected=0.005, m_isolated=0.0005,
                cycle_len=100_000, warm_len=20_000, scenario="3peak")
```

---

## Output

Each replicate is one `.rds` holding the list returned by `run_both()`:

```
$warm , $cold                 one entry per thermal ecology
   $demog                     climate, occupancy and per-deme Ne through time
   $gl_present                genlight, 20 diploids per occupied deme at t = 0
   $pi                        data frame, one row per sampling time:
                                time, pi, tajima_d, theta_w, fst_mean, he_mean
   $gl_timeseries             genotypes at every sampling time (NULL by default)
$meta                         scenario, replicate, seed and all parameters
```

Alongside them, `design_matrix.csv` records what was run and `session_info.txt`
records the R and Python versions that ran it.

### Two different Fst values, deliberately

They are not interchangeable, and they answer different questions:

- **`$pi$fst_mean`** — computed inside the simulation at every sampling time,
  averaged over **four-neighbour (adjacent) deme pairs only**. This is the
  local-structure time series shown in Figure 2B.
- **`pw$fst`** — recomputed by `R/05_extract_fst.R` from the present-day
  genlight objects, over **all pairs of demes**, with the geographic distance
  between them. Needed for Figures 5 and 6, where isolation by distance has to
  be separated from the barrier effect.

Both use Hudson's estimator as a ratio of averages, the low-bias choice when
demes differ in sample size and most SNPs are near-monomorphic.

---

## Grid coordinates

Populations are named `p<k>` with `k = row * n + column`, both 0-based, matching
`cell_idx()` in the Python engine. `R/05_extract_fst.R` converts these to
1-based `row` and `col` columns. The river sits on Python row `(n %/% 2) - 1`,
which is `n %/% 2` in those 1-based tables — that is how
`R/07_fig_h2_h3.R` decides whether a pair of demes spans it. Pairs with a deme
sitting *on* the river row cannot be assigned to a bank and are dropped rather
than forced to one side.

---

## An optional dependency

`plot_fst_present_scenarios()` in `R/09_fig_supplementary.R` is the only
function in the repository that needs `dartRverse` (it calls
`dartR.base::gl.report.fstat()`). Nothing in the main pipeline does: the
present-day Fst behind Figures 5 and 6 is computed independently in
`R/05_extract_fst.R`. You can reproduce every published figure without
installing `dartRverse`.

## Citation

msprime: Baumdicker et al. (2022) *Genetics* 220: iyab229.
