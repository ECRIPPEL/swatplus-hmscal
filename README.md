# swatplus-hmscal

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19793606.svg)](https://doi.org/10.5281/zenodo.19793606)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**HMSCal: Hierarchical Multiscale Calibration of SWAT+ guided by Adaptive Sensitivity Analysis.**

A reproducible R-based framework for calibrating SWAT+ hydrological models across three temporal scales (annual water balance → monthly streamflow → daily streamflow), using Morris sensitivity screening re-executed at each iteration, σ-informed parameter-range narrowing, and cross-iteration filtering to manage equifinality. Built on top of [SWATrunR](https://chrisschuerz.github.io/SWATrunR/).

What separates HMSCal from a single-scale workflow is that a runoff coefficient target and a multisite aquifer BFI proxy stay active in **every** iteration, not only at the end. A parameter set advances only when it meets every criterion of the current iteration, so the accepted ensemble is constrained on internal flow partitioning not only on outlet streamflow.

If you use this software in research, please see [Citation](#citation) below.

> **Note on the name.** This project was released as `swatplus-hms-calibration` up to v1.0.0. From v1.1.0 it is named **HMSCal** (`swatplus-hmscal`), to avoid collision with HEC-HMS and to make clear that the acronym denotes a calibration strategy, not a hydrological model. The old repository URL still redirects, and the v1.0.0 tag and its Zenodo deposit are unchanged.

---

## ⚠️ Methodology note — the HMSCal acceptance criterion (RC + multisite BFI)

The framework uses a **joint acceptance criterion** that stays active through all four iterations:

```
satisfactory_run = (obj_rc      <= threshold_obj_rc)         AND
                   (obj_bfi_mdae <= threshold_obj_bfi)         AND
                   (NSE          >= threshold_nse,    Iterations 2 to 4 only) AND
                   (KGE          >= threshold_kge,    Iterations 2 to 4 only) AND
                   (|PBIAS|      <= threshold_abs_pbias, Iterations 2 to 4 only)
```

where:
- `obj_rc       = |RC_obs − RC_sim|`, with `RC = Q_total / P` (single dimensionless scalar)
- `obj_bfi_mdae = median over N gauges of |BFI_obs_i − BFI_sim_i|`

The two sides of that comparison are deliberately **not** the same estimator, and this is the point of the method.

- `BFI_obs_i` is the **reference slow-flow signature**, obtained by applying Eckhardt's (2005) forward digital filter to the *observed* daily series at gauge `i`, with `(a, BFImax)` from the Vogel & Kroll recession analysis and the Collischonn & Fan q90/q50 regression.
- `BFI_sim_i` is the **SWAT+ internal aquifer partition**: cumulative aquifer-to-channel discharge over the upstream floodplain aquifers, divided by cumulative channel discharge at that gauge. It is read from the model's internal fluxes, not from the shape of the simulated hydrograph.

Comparing them assumes the reference slow flow should be carried mainly by the SWAT+ aquifer pathway. Filtering on the internal partition is what forces calibration to actually activate the aquifer module, instead of mimicking a baseflow recession with slow lateral flow. A digital filter applied to *simulated* flow cannot tell those two apart, because it sees only hydrograph shape.

This is set by `bfi.sim_source: "swat_internal_bf"` in every config. The alternative, `sim_source: "eckhardt"`, applies the filter to simulated flow and is kept as a **diagnostic only**; it is not the acceptance criterion. See [R/bfi_eckhardt.R](R/bfi_eckhardt.R).

The proxy error is the median absolute deviation across the gauges that have a directly mapped aquifer discharge. In the published application that is three of the nine gauges; the remaining six are used to check spatial consistency, not as calibration targets.


### Required outputs in `cal_outputs`

To compute `RC_sim`:
```yaml
cal_outputs:
  wateryld:  { file: "basin_wb_aa",  variable: "wateryld", unit: 1 }   # surq + latq
  baseflow:  { file: "basin_aqu_aa", variable: "flo_cha",  unit: 1 }   # gw_q
  precip:    { file: "basin_wb_aa",  variable: "precip",   unit: 1 }
  et:        { file: "basin_wb_aa",  variable: "et",       unit: 1 }   # diagnostic only
```

To compute `BFI_sim_i`:
```yaml
cal_outputs:
  vazao_multi_day:                  # daily Q at all selected gauges
    file: "channel_sdmorph_day"
    variable: "flo_out"
    unit: [5, 19, 42, 73, 79, 114, 136, 146, 194]
```

### Required reference table

The CSV with per-gauge `BFI_obs`, recession constant `a`, and Eckhardt `BFImax` must be placed at `<project_path>/<bfi_targets_csv>`. Default name `bfi_targets.csv`. Required columns:

| column         | meaning                                         |
|----------------|-------------------------------------------------|
| `channel_id`   | SWAT+ unit ID of the channel                    |
| `bfi`          | observed BFI (target, dimensionless [0,1])      |
| `a`            | recession constant from observed Q (e.g., q90/q50) |
| `bfimax_q9050` | Eckhardt asymptotic BFI (Collischonn & Fan)     |

A template is shipped at [data/bfi_targets.csv](data/bfi_targets.csv) — copy it to your project_path and edit per your gauges.

---

## Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Project Structure](#2-project-structure)
3. [Before You Start](#3-before-you-start)
4. [Understanding the Workflow](#4-understanding-the-workflow)
5. [Configuration Files (YAML)](#5-configuration-files-yaml)
6. [Running the Full Pipeline](#6-running-the-full-pipeline)
7. [Running Individual Steps](#7-running-individual-steps)
8. [Resuming After Interruption](#8-resuming-after-interruption)
9. [Understanding the Outputs](#9-understanding-the-outputs)
10. [Adapting to Your Catchment](#10-adapting-to-your-catchment)
11. [Troubleshooting](#11-troubleshooting)

[Release history](#release-history) · [Citation](#citation) · [Contact](#contact) · [License](#license)

---

## 1. Prerequisites

### Software

| Software   | Version  | Notes                                                                        |
|------------|----------|------------------------------------------------------------------------------|
| R          | >= 4.3.3                                 | Required for SWATrunR compatibility          |
| SWAT+      | rev. 2026.61.0.2.61 (`swatplus-61.0.2.61-ifx-win_amd64-Rel`) | The revision this was built and tested against. Executable must be inside project_path. Other revisions may write output files with a different shape, see [`SWATrunR_mod/`](SWATrunR_mod/) |
| RTools     | >= 4.3 (optional)                        | Windows only. Not needed for SWATrunR, which is pure R; only if you build some other package from source |

### R Packages

Install all required packages before running:

```r
install.packages(c(
  "yaml",         # YAML config parser
  "dplyr",        # Data manipulation
  "tibble",       # Modern data frames
  "tidyr",        # Data reshaping
  "purrr",        # Functional programming
  "hydroGOF",     # NSE, KGE, PBIAS functions
  "lubridate",    # Date handling
  "ggplot2",      # Plotting
  "ggrepel",      # Label placement in plots
  "lhs"           #Latin Hypercube Sampling
))
```

#### SWATrunR needs a patched install

`SWATrunR` is not on CRAN, and a plain `remotes::install_github("chrisschuerz/SWATrunR")` is **not enough** to run this pipeline. Stock SWATrunR 0.9.4 hits two problems: it fails to install on readr 2.2.0 and newer, and it aborts `run_swatplus()` after every SWAT+ thread has already finished whenever an `lsunit_wb_aa` output is requested, destroying the results.

Both problems, their causes, and the patch are documented in **[`SWATrunR_mod/`](SWATrunR_mod/)**. Go through that folder before your first run.

From R, with this repository's root as your working directory:

```r
source("SWATrunR_mod/install_swatrunr_patched.R")
install_swatrunr_patched()
```

The script pins `readr`, downloads SWATrunR at the pinned upstream commit, applies the fix, installs it, and verifies the installed package before reporting success. Restart R afterwards. No git and no compiler needed: SWATrunR is a pure R package. See [`SWATrunR_mod/README.md`](SWATrunR_mod/README.md) for the manual route and for what the fix actually does.

### SWAT+ Project

Your SWAT+ TxtInOut folder (referred to as `project_path` in configs) must contain:
- A working SWAT+ model that runs successfully
- `cal_parms.cal` (used for parameter range validation in Iteration 1)
- `soils.sol`, `hydrology.hyd`, `topography.hyd` (read during range checks)

### Observed Data

You need CSV files with observed streamflow placed inside your `project_path`:

| File                   | Columns       | Format              | Used in        |
|------------------------|---------------|---------------------|----------------|
| `obs_flow_monthly.csv` | Date, Flow    | YYYY-MM-DD, m^3/s   | Iterations 1 and 2    |
| `obs_flow_daily.csv`   | Date, Flow    | YYYY-MM-DD, m^3/s   | Iteration 3        |

- `Date` must be in `YYYY-MM-DD` format (e.g., `1980-01-01`)
- `Flow` is streamflow in m^3/s
- Monthly file: one row per month (use the 1st day of each month as date)
- Daily file: one row per day

#### Example data shipped with the repository

[`data/obs/`](data/obs) carries a worked example from the Upper Uruguay application: the basin outlet plus the three gauges that define the multisite aquifer partition proxy, daily and monthly.

| Channel | ANA gauge | Role in the example                        |
|---------|-----------|--------------------------------------------|
| 33      | 74100000  | Iraí, basin outlet (streamflow calibration target) |
| 79      | 73350000  | BFI reference gauge                        |
| 194     | 70200000  | BFI reference gauge                        |
| 73      | 71498000  | BFI reference gauge                        |

[`data/obs/obs_index.csv`](data/obs/obs_index.csv) records the metadata for these four (coordinates, quality class, period of record, regulation status). The three BFI gauges are the ones listed in [`data/bfi_targets.csv`](data/bfi_targets.csv).

This is an **example subset**. The complete series come from Brazil's National Water Agency (ANA) Hidroweb portal, <https://www.snirh.gov.br/hidroweb/>. Point `obs_flow_daily.csv` / `obs_flow_monthly.csv` at your own gauges to calibrate another catchment.

---

## 2. Project Structure

```
swatplus-hmscal/
├── config/
│   ├── iteration1.yaml     # Iteration 1: annual water balance      (MORRIS1 + ITER1)
│   ├── iteration2.yaml     # Iteration 2: monthly streamflow        (MORRIS2 + ITER2)
│   ├── iteration3.yaml     # Iteration 3: daily streamflow          (MORRIS3 + ITER3)
│   ├── iteration4.yaml     # Iteration 4: daily, narrowed ranges    (ITER4, no new screening)
│   └── validation.yaml     # Split-sample validation on an independent decade
├── R/
│   ├── baseflow_comparison.R # Observed vs simulated baseflow diagnostics
│   ├── bfi_eckhardt.R      # Eckhardt (2005) digital filter + Collischonn & Fan BFImax
│   ├── build_outputs.R     # Assembles the cal_outputs request sent to SWATrunR
│   ├── check_range.R       # Parameter range validation vs cal_parms.cal
│   ├── lhoat_engine.R      # LH-OAT trajectory generation + Elementary Effects
│   ├── metrics.R           # NSE, KGE, PBIAS, RC, BFI and water balance metrics
│   ├── morris_classify.R   # Morris screening classification + range narrowing
│   ├── param_inherit.R     # Cross-iteration parameter inheritance (drop / override / extras)
│   ├── plots.R             # Publication-ready scientific plots
│   ├── run_filter.R        # Invalid run detection and filtering
│   ├── versioning.R        # Auto-versioned output folder management
│   └── write_run_info.R    # Generates the run_info.txt metadata file
├── data/
│   ├── bfi_targets.csv       # Per-gauge BFI_obs, recession constant a, BFImax
│   ├── bfi_targets_1990s.csv # Same, re-derived for the validation decade
│   └── obs/                  # Example observed streamflow (see Observed Data)
├── SWATrunR_mod/           # READ FIRST: two SWATrunR defects that block this
│   ├── README.md           #   pipeline, their causes, and how to install
│   ├── install_swatrunr_patched.R          # one-command patched install
│   └── 01-read_output_i-header-guard.patch # the diff, for manual use
├── master.R                # Orchestrator: runs all 7 steps sequentially
├── run_morris.R               # Generic Morris screening script (works for any iteration)
├── run_cal.R               # Generic calibration script (works for any iteration)
└── run_val.R               # Split-sample validation of an accepted ensemble
```

**You should NOT need to edit any `.R` file.** All customization is done through the YAML config files.

---

## 3. Before You Start

### Set Your Working Directory

All scripts expect to run from the project root (`swatplus-hmscal/`).

**In RStudio:**
```r
setwd("C:/path/to/your/swatplus/txtinout/swatplus-hmscal")
```

**In VSCode (R terminal):**
```r
setwd("C:/path/to/your/swatplus/txtinout/swatplus-hmscal")
```

**From command line:**
```bash
cd C:/path/to/your/swatplus/txtinout/swatplus-hmscal
```

### Verify Your SWAT+ Model Runs

Before calibration, ensure your SWAT+ model runs successfully:

```r
library(SWATrunR)
sim_test <- run_swatplus(
  project_path = "C:/path/to/your/swatplus/txtinout",
  output = define_output(file = "channel_sd", variable = "flo_out", unit = 33),
  start_date = "1979-01-01",
  end_date   = "1980-12-31",
  years_skip = 1,
  n_thread   = 4
)
```

If this fails, fix your SWAT+ model before proceeding.

---

## 4. Understanding the Workflow

HMSCal runs **four iterations** that progress from annual water-balance constraints to monthly and then daily streamflow criteria. A Morris screening precedes Iterations 1 to 3; Iteration 4 reuses the classification obtained from the third screening. That makes 7 executable steps:

```
MORRIS1  ->  ITER1  ->  MORRIS2  ->  ITER2  ->  MORRIS3  ->  ITER3  ->  ITER4
(step 1)   (step 2)   (step 3)    (step 4)   (step 5)    (step 6)   (step 7)
```

Each iteration draws 250 parameter sets by Monte Carlo sampling from the current ranges. A simulation is **accepted** only when it satisfies *every* criterion active in that iteration. The accepted sets, together with the Morris classification, define the ranges sampled in the next iteration. No parameter is adjusted by hand between iterations.

Each step writes to its own auto-versioned folder named `<TAG>_v<N>/` (e.g. `MORRIS1_v1/`, `ITER4_v3/`). All artefacts inside follow the pattern `<descriptor>_<TAG>.<ext>`, so any file stays self-identifying when copied out.

**What stays active throughout.** The runoff coefficient (RC) criterion and the multisite BFI proxy criterion are evaluated in *all four* iterations. They are not an opening stage that later gets dropped: that persistence is what constrains the internal flow partition rather than the outlet hydrograph alone. What changes between iterations is the streamflow criterion layered on top, and the temporal scale it is evaluated at.

### Iteration 1: annual water balance
- **Tags:** MORRIS1 (screening), ITER1 (calibration)
- **Scale:** annual
- **Criteria:** `|RC_ref - RC_sim| <= 0.10` and `BFI proxy MdAE <= 0.10`
- **Parameters:** 14 (*esco*, *epco*, *perco*, *latq_co*, *cn3_swf*, *revap_co*, *canmx*, *petco*, *bf_max*, *flo_min*, *revap_min*, *alpha_bf*, *cn2*, *awc*)
- **Screening:** 15 LH-OAT Morris trajectories, `r(p+1)` = 225 runs
- **Output:** `inherited_ranges_ITER1_v<N>.csv`, consumed by MORRIS2

### Iteration 2: monthly streamflow
- **Tags:** MORRIS2, ITER2
- **Scale:** monthly
- **Criteria:** monthly NSE >= 0.70, KGE >= 0.70, |PBIAS| <= 10%, **plus** the RC and BFI proxy criteria
- **Parameters:** the same 14, on the ranges narrowed after Iteration 1
- **Screening:** 225 runs
- **Output:** `inherited_ranges_ITER2_v<N>.csv`, consumed by MORRIS3

### Iteration 3: daily streamflow
- **Tags:** MORRIS3, ITER3
- **Scale:** daily
- **Criteria:** daily NSE >= 0.50, KGE >= 0.50, |PBIAS| <= 15%, **plus** the RC and BFI proxy criteria
- **Parameters:** 18. The four routing parameters (*surlag*, *lat_ttime*, *chn*, *ovn*) enter here, because they only matter once the criterion is evaluated daily
- **Screening:** `r(p+1)` = 285 runs
- **Output:** `inherited_ranges_ITER3_v<N>.csv`, consumed by ITER4

### Iteration 4: daily criteria on the narrowed ranges
- **Tags:** ITER4 (no new screening)
- **Scale:** daily
- **Criteria:** daily NSE >= 0.50, KGE >= 0.50, |PBIAS| <= 15%, **plus** the RC and BFI proxy criteria
- **Parameters:** the same 18, on the ranges narrowed after Iteration 3
- **Screening:** none. Iteration 4 reuses the MORRIS3 classification, which is why `config/iteration4.yaml` still points `morris.tag_prefix` at `MORRIS3`
- **Output:** the final accepted ensemble, and the input to `run_val.R`

### How the iterations connect

Each iteration inherits parameter ranges from the previous one through auto-discovered files:

1. **ITER1** writes `ITER1_v<N>/inherited_ranges_ITER1_v<N>.csv`. **MORRIS2** and **ITER2** locate the latest `ITER1_v*` folder via `inherit.from: "ITER1"`.
2. **ITER2** writes `inherited_ranges_ITER2_v<N>.csv`, consumed via `inherit.from: "ITER2"`.
3. **ITER3** writes `inherited_ranges_ITER3_v<N>.csv`, consumed by **ITER4** via `inherit.from: "ITER3"`.

Filenames are fully determined by the folder TAG. There is no `inherit.ranges_file` field.

> **A note on two names in the code.** The YAML key `iteration:` is the iteration number, 1 to 4. Separately, some generated artefacts and one config value still use the word *satisfactory* (`cal.inheritance_source: "satisfactory"`, the `satisfactory_results` field, `simulation_results_satisfactories_<TAG>.csv`). Those are literal identifiers the code reads and writes, so they are left untouched. In this document and in the paper, a set that clears every active criterion is called **accepted**.

### Inheritance source — how the CSV is filled

The `inherited_ranges_<TAG>.csv` is computed via Morris-stratified quantile narrowing (`calc_new_ranges()`). The set of runs that feeds it is controlled by `cal.inheritance_source` in the YAML:

- `"satisfactory"` (default) — only accepted runs feed the narrowing. If no run is accepted, the producer **stops loud** with an actionable error message instead of silently skipping the file.
- `"all_valid"` — all valid runs feed the narrowing (the resulting ranges are not filtered on acceptance, so this is exploratory mode for keeping the pipeline moving).

Parameter ranges are narrowed using quantile rules based on Morris sensitivity groups:
- **High importance + interaction** (mu* > mean, sigma > mean): Q5–Q95 (wide exploration)
- **High importance** (mu* > mean, sigma <= mean): Q10–Q90 (moderate narrowing)
- **Low importance** (mu* <= mean): Q25–Q75 (aggressive narrowing)

### Inheritance flexibility — drop / override / extras

The `inherit:` block supports three optional operations applied in order:

```yaml
inherit:
  from: "ITER1"

  # (a) Drop parameters that became uninformative — fix at default
  drop: ["revap_min", "flo_min"]

  # (b) Override min/max of inherited (or extra) parameters
  override_ranges:
    cn2:   { min: -40, max: 25 }
    petco: { min: 0.85, max: 1.15 }

extra_params:                          # (c) Add new parameters not present in the source iteration
  - short: "surlag"
    par_name: "surlag::surlag.bsn | change = absval"
    min: 1
    max: 24
```

All three are validated at config-read time — referencing a parameter that doesn't exist in the inherited set or in `extra_params` produces an immediate `stop()` with a clear message.

---

## 5. Configuration Files (YAML)

Each YAML config controls one iteration. Here are the key sections you may want to customize:

### General Settings

```yaml
iteration: 1                                              # Iteration number (1 to 4)
project_path: "C:/path/to/your/swatplus/txtinout"     # Path to your SWAT+ TxtInOut folder
seed: 123                                             # Random seed for reproducibility
threads: 16                                           # Number of parallel SWAT+ threads
script_version: "1.0.0"                               # Version tag saved with results
```

### Morris screening settings

```yaml
morris:
  tag_prefix: "MORRIS1"                # Output folder pattern: MORRIS1_v<N> (auto-versioned)
  m_traj: 1                         # Number of LH-OAT trajectories
  delta: 0.20                       # Perturbation step size (fraction of range)
  start_date: "1979-01-01"          # screening simulation period
  end_date: "1980-12-31"
  warmup: 1                         # Warmup years (skipped in analysis)
```

The `tag_prefix` values follow the master.R step naming: `MORRIS1`, `MORRIS2`, `MORRIS3` for the three iterations respectively.

> **Tip:** For Iteration 1 screening, `m_traj: 1` is often enough. For Iterations 2 to 4, use
> `m_traj: 15` or more for robust Morris indices. More trajectories = more SWAT+ runs
> but more reliable sensitivity ranking. Total runs = m_traj × (p + 1), where p is
> the number of parameters.

### Calibration Settings

**`type: "soft"` (Iteration 1):**
```yaml
cal:
  type: "soft"                      # RC + BFI proxy criteria (Iteration 1)
  tag_prefix: "ITER1"                # Output folder pattern: ITER1_v<N>
  n_simulations: 250                # Number of random parameter sets
  inheritance_source: "satisfactory"  # "satisfactory" (default, fail loud if 0) | "all_valid"
  threshold_pct: 15                 # Maximum relative error (%) for water balance
  targets:
    et_rto: 0.48147                 # Target ET/P ratio
    wyld_rto: 0.51853               # Target WYLD/P ratio (= 1 - ET/P)
```

**`type: "hard"` (Iterations 2 to 4):**
```yaml
cal:
  type: "hard"                      # Performance thresholds
  tag_prefix: "ITER2"                # or "ITER3", "ITER4"
  n_simulations: 250
  inheritance_source: "satisfactory"
  threshold_nse: 0.70               # Minimum NSE
  threshold_kge: 0.70               # Minimum KGE
  threshold_abs_pbias: 10           # Maximum |PBIAS| (%)
```

**Diagnostic-only WB closure (the daily iterations)** — declared in `iteration3.yaml` and `iteration4.yaml`:
```yaml
diagnostic_balance:
  et_rto: 0.48147               # target ET/P (annual closure)
  wyld_rto: 0.51853             # target Q_total/P
  threshold_pct: 15             # legacy field; unused in the daily iterations (diagnostic only)
```
This block triggers `calc_balance_diagnostic()` to attach `rto_et`, `rto_q`, `rto_wyld`, `rto_baseflow` columns to every daily-iteration run + a `balance_scatter_<TAG>.tif` figure. **It does NOT contribute to `satisfactory_run`** — that is defined solely by daily NSE/KGE/PBIAS.

The `cal.inheritance_source` flag drives which set of runs feeds the narrowing that produces `inherited_ranges_<TAG>.csv`:

| Value | Source | If empty |
|---|---|---|
| `"satisfactory"` (default) | `simulation_results_satisfactories_<TAG>.csv` | `stop()` LOUD with actionable message |
| `"all_valid"` | `simulation_results_<TAG>.csv` (all valid runs) | only fails if 0 valid runs total (rare) |

### Parameter definitions (Iteration 1 only)

Iteration 1 defines parameters explicitly. Each parameter needs:

```yaml
parameters:
  - short: "cn2"                                        # Short name (used in plots)
    par_name: "cn2::cn2.hru | change = pctchg"          # SWATrunR notation
    min: -15                                             # Lower bound
    max: 15                                              # Upper bound
```

**SWATrunR parameter notation:**
- `variable::file | change = type`
- Change types: `absval` (absolute value), `pctchg` (% change), `abschg` (absolute change)
- Example: `awc::soils.sol | change = pctchg` means modify AWC in soils.sol by a percentage

### Parameter inheritance (Iterations 2 to 4)

Iterations 2 to 4 inherit ranges from the previous iteration via TAG-based auto-discovery — the file name is fully determined by the source folder TAG, so no `ranges_file` field is needed:

```yaml
inherit:
  from: "ITER1"                      # tag prefix of source iteration (no _vN)

  # Optional: drop parameters from inherited set (fix at default)
  # drop: ["revap_min", "flo_min"]

  # Optional: override min/max of an inherited (or extra) parameter
  # override_ranges:
  #   cn2:   { min: -40, max: 25 }
  #   petco: { min: 0.85, max: 1.15 }

param_map:
  cn2:    "cn2::cn2.hru | change = pctchg"
  awc:    "awc::soils.sol | change = pctchg"
  # ... maps short names to SWATrunR notation
```

Source TAGs by iteration:
- `iteration2.yaml` → `inherit.from: "ITER1"`
- `iteration3.yaml` → `inherit.from: "ITER2"`
- `iteration4.yaml` → `inherit.from: "ITER3"`

Each pulls the latest `<TAG>_v<N>/inherited_ranges_<TAG>_v<N>.csv` from `project_path`.

### Adding extra parameters in later iterations

Iteration 3 introduces new parameters not present in Iteration 2:

```yaml
extra_params:
  - short: "surlag"
    par_name: "surlag::surlag.bsn | change = absval"
    min: 1
    max: 24
```

`extra_params` is also subject to `inherit.override_ranges` — you can add a parameter and immediately reposition its prior in the same yaml.

### Extra Outputs (per-channel / per-LSU / per-HRU / per-aquifer)

The base `cal_outputs:` / `morris_outputs:` blocks enumerate one entry per (file, variable, unit) — practical for a handful of basin-aggregated outputs but unwieldy for per-channel, per-LSU, per-HRU or per-aquifer retention (hundreds of units).

For bulk retention, use the optional `extra_outputs:` block. It mirrors the `h3_outputs:` pattern from the sister project [swatplus-calibrator](https://github.com/swat-model/swatplus-calibrator) and emits one SWATrunR `define_output()` per variable with a vector of units behind the scenes:

```yaml
extra_outputs:
  enabled: true                  # master switch

  # channel_sd_aa per reach — for channel dissipation (H_diss) + per-channel BFI
  channels:
    enabled: true
    units: [22, 33, 79, 113, 209]   # OR units: "auto" + cha_units_file
    vars: ["flo_in", "flo_out", "evap", "seep"]

  # lsunit_wb_aa per LSU — spatial water balance decomposition
  lsu:
    enabled: false
    n_lsu: 422
    wb_vars: ["precip", "et", "wateryld", "surq_gen", "latq", "perc"]

  # hru_wb_aa per filtered HRU — per-class LULC banding
  hru:
    enabled: false
    hru_ids_file: "data/hru_ids_for_h3.txt"   # one int per line
    wb_vars: ["precip", "et", "surq_gen", "latq", "perc", "wateryld"]

  # aquifer_aa per aquifer — per-aquifer water balance + baseflow
  aquifer:
    enabled: false
    n_aqu: 47
    aqu_vars: ["flo", "dep_wt", "stor", "rchrg", "seep", "revap",
               "flo_cha", "flo_res", "flo_ls"]
```

Each sub-block (`channels`, `lsu`, `hru`, `aquifer`) is independently toggleable via its own `enabled:` field. The master `extra_outputs.enabled: true` is required for any sub-block to be processed.

**Where outputs land** — once captured, the variables are accessible in the saved RDS as additional entries inside `sim$simulation`:

| extra_outputs block | Naming pattern in `sim$simulation`             |
|---------------------|------------------------------------------------|
| `channels`          | `cha_aa_<var>` (e.g., `cha_aa_flo_out`)        |
| `lsu`               | `lsu_<var>` (e.g., `lsu_precip`)               |
| `hru`               | `hru_<var>` (e.g., `hru_et`)                   |
| `aquifer`           | `aqu_<var>` (e.g., `aqu_rchrg`)                |

Each entry is a tibble with `date` plus one column per unit (e.g., `flo_out_22`, `flo_out_33`, …).

**Memory cost** — roughly 10-50 MB per simulation depending on which sub-blocks are enabled. Practical guidance:

- **Morris screening steps (high run count):** keep `extra_outputs.enabled: false` to keep the per-thread SQLite small.
- **Diagnostic CAL iterations (e.g., ITER4 BFI multi-site, channel-dissipation analysis):** enable the specific sub-blocks you need.
- **Full retention pass:** to mirror the calibrator's `validation_multi`/`h3_outputs` pattern (run the accepted ensemble once with everything captured), enable all four sub-blocks in a dedicated calibration YAML.

This block is OPTIONAL — if absent or `enabled: false`, behaviour is identical to the legacy flat `cal_outputs:` / `morris_outputs:` schema.

---

## 6. Running the Full Pipeline

### Option A: Run Everything from Step 1

```r
setwd("C:/path/to/your/swatplus/txtinout/swatplus-hmscal")
source("master.R")
```

This will execute all 7 steps sequentially:
```
Step 1/7 - MORRIS1 - Morris screening, annual scale
Step 2/7 - ITER1 - Iteration 1, annual water balance (RC + BFI proxy)
Step 3/7 - MORRIS2 - Morris screening, monthly scale
Step 4/7 - ITER2 - Iteration 2, monthly streamflow + RC + BFI proxy
Step 5/7 - MORRIS3 - Morris screening, daily scale
Step 6/7 - ITER3 - Iteration 3, daily streamflow + RC + BFI proxy (routing params enter)
Step 7/7 - ITER4 - Iteration 4, same daily criteria on the narrowed ranges
```

Each step prints timing information when completed.

### Option B: From the Command Line

```bash
cd C:/path/to/your/swatplus/txtinout/swatplus-hmscal
Rscript master.R
```

---

## 7. Running Individual Steps

You can run any step independently.

### Morris screening (any iteration)

**From R console / RStudio / VSCode:**
```r
ITERATION_CONFIG <- "config/iteration1.yaml"
source("run_morris.R")
```

**From command line:**
```bash
Rscript run_morris.R config/iteration1.yaml
Rscript run_morris.R config/iteration2.yaml
Rscript run_morris.R config/iteration3.yaml
```

### Calibration (any iteration)

**From R console / RStudio / VSCode:**
```r
ITERATION_CONFIG <- "config/iteration2.yaml"
source("run_cal.R")
```

**From command line:**
```bash
Rscript run_cal.R config/iteration1.yaml
Rscript run_cal.R config/iteration2.yaml
Rscript run_cal.R config/iteration3.yaml
```

> **Important:** Calibration steps depend on the Morris screening of the same iteration.
> Always run the screening before the calibration of each iteration. Iteration 2+ also depends on calibration
> results from the previous iteration (for parameter range inheritance).

### Correct Execution Order

If running steps individually, follow this order:

```
1. Rscript run_morris.R  config/iteration1.yaml        # MORRIS1
2. Rscript run_cal.R  config/iteration1.yaml        # ITER1   (needs MORRIS1)
3. Rscript run_morris.R  config/iteration2.yaml        # MORRIS2   (needs ITER1 inherited_ranges)
4. Rscript run_cal.R  config/iteration2.yaml        # ITER2   (needs MORRIS2)
5. Rscript run_morris.R  config/iteration3.yaml        # MORRIS3   (needs ITER2 inherited_ranges)
6. Rscript run_cal.R  config/iteration3.yaml        # ITER3 (needs MORRIS3)
7. Rscript run_cal.R  config/iteration4.yaml  # ITER4 (needs ITER3 inherited_ranges)
```

---

## 8. Resuming After Interruption

If the pipeline is interrupted (crash, power loss, manual stop), you can resume from any step:

```r
START_FROM <- 3        # Resume from the Iteration 2 screening (step 3)
source("master.R")
```

Steps before `START_FROM` will be skipped with `[SKIP]` in the console output.

### Step Reference

| Step | Tag | Script | Config | What it does |
|---|---|---|---|---|
| 1 | MORRIS1 | `run_morris.R` | `config/iteration1.yaml` | Morris screening at the annual scale |
| 2 | ITER1 | `run_cal.R` | `config/iteration1.yaml` | **Iteration 1**: RC and BFI proxy criteria |
| 3 | MORRIS2 | `run_morris.R` | `config/iteration2.yaml` | Morris screening at the monthly scale |
| 4 | ITER2 | `run_cal.R` | `config/iteration2.yaml` | **Iteration 2**: monthly NSE/KGE/PBIAS, plus RC and BFI proxy |
| 5 | MORRIS3 | `run_morris.R` | `config/iteration3.yaml` | Morris screening at the daily scale |
| 6 | ITER3 | `run_cal.R` | `config/iteration3.yaml` | **Iteration 3**: daily NSE/KGE/PBIAS, plus RC and BFI proxy. Routing parameters enter |
| 7 | ITER4 | `run_cal.R` | `config/iteration4.yaml` | **Iteration 4**: same daily criteria on the ranges narrowed after Iteration 3. Reuses the MORRIS3 classification |

Validation is run separately, after the four iterations: `run_val.R` with `config/validation.yaml` propagates the Iteration 4 accepted ensemble through an independent decade.

> **Tip:** If a step fails, fix the issue and resume from that step. You don't need
> to re-run previous steps — their results are already saved in versioned folders.

---

## 9. Understanding the Outputs

### Output Folder Structure

Each step creates an auto-versioned folder inside your `project_path`. The folder name follows the pattern `<TAG>_v<N>/`, and every artefact inside follows the pattern `<descriptor>_<TAG>.<ext>`:

```
C:/path/to/your/swatplus/txtinout/
├── MORRIS1_v1/               # MORRIS1 results (1st run)
├── MORRIS1_v2/               # MORRIS1 results (2nd run, if re-run)
├── ITER1_v1/               # ITER1 (soft) results
├── MORRIS2_v1/               # MORRIS2 results
├── ITER2_v1/               # ITER2 (hard, monthly) results
├── MORRIS3_v1/               # MORRIS3 results
├── ITER3_v1/             # ITER3 (hard, daily Iter 1)
└── ITER4_v1/             # ITER4 (hard, daily Iter 2 — depth-refined ensemble)
```

Version numbers auto-increment. If you re-run a step, a new folder (v2, v3...) is created — previous results are never overwritten.

### Morris screening output files (inside each `MORRIS<N>_v<V>/`)

| File                                  | Description                                                |
|---------------------------------------|------------------------------------------------------------|
| `run_info.txt`                        | Human-readable metadata (period, params, targets, counts)  |
| `morris_results_<TAG>.rds`            | Full R object with all Morris screening results                         |
| `morris_ranking_<TAG>.csv`            | Morris sensitivity ranking (mu*, sigma per param)          |
| `morris_metrics_<TAG>.csv`            | Performance metrics for all LH-OAT runs                    |
| `morris_ee_<TAG>.csv`                 | Elementary Effects for all trajectories                    |
| `param_info_<TAG>.csv`                | Parameter info (names, ranges) used in the screening                 |
| `morris_screening_<TAG>.tif`          | Morris plot (mu* vs sigma) — publication-ready             |
| `check_param_range_<TAG>.csv`         | Parameter range validation (MORRIS1 only)                     |

`<TAG>` is the folder name (e.g., `MORRIS1_v3`, `MORRIS3_v1`). Filenames are deterministic so a CSV can be located unambiguously by tag.

### Iteration 4 outputs (depth refinement)

Iterations 3 and 4 run at the same daily scale. Iteration 4 is a depth refinement: it re-samples from the ranges narrowed by the sets Iteration 3 accepted, and reuses the MORRIS3 classification instead of screening again.

1. **ITER3** (step 6) — `run_cal.R` + `config/iteration3.yaml`. Daily NSE/KGE/PBIAS, plus the RC and BFI proxy criteria. Writes `ITER3_v<N>/inherited_ranges_ITER3_v<N>.csv` from the accepted subset (Morris-stratified quantile narrowing).

2. **ITER4** (step 7) — `run_cal.R` + `config/iteration4.yaml`. `inherit.from: "ITER3"`. Auto-discovers `ITER3_v<N>/inherited_ranges_ITER3_v<N>.csv` and re-samples a fresh LHS from those narrower ranges. Same daily criteria, deeper sampling around the accepted region.

**Multi-scale separation premise.** Both iterations calibrate ONLY on daily NSE/KGE/PBIAS. The WB closure quantities (rto_et, rto_q, rto_wyld, rto_baseflow) are computed and reported as DIAGNOSTIC columns + figure on every daily-iteration run, but never combined with the daily filter — combining them would mix scales and undermine the multi-scale separation that defines this framework.

When run via `source("master.R")`, both steps execute automatically. To run them standalone:

```r
# Iteration 1 (step 6)
ITERATION_CONFIG <- "config/iteration3.yaml"
source("run_cal.R")
# → produces ITER3_v<N>/inherited_ranges_ITER3_v<N>.csv

# Iteration 2 (step 7) — uses iteration4.yaml
ITERATION_CONFIG <- "config/iteration4.yaml"
source("run_cal.R")
# → produces ITER4_v<K>/  (depth-refined daily ensemble)
```

Both iterations require `baseflow` in `cal_outputs` (basin_aqu_aa.flo_cha) so the diagnostic `Q_total = wateryld + baseflow` formulation can be reported. See the methodology note at the top of this README.

### Calibration output files (inside each `ITER<N>_v<V>/`)

| File                                                | Description                                                              |
|-----------------------------------------------------|--------------------------------------------------------------------------|
| `run_info.txt`                                      | Human-readable metadata (period, params, targets, counts, inherit info)  |
| `simulation_results_<TAG>.rds`                      | Full R object: config, sim, results, satisfactory_results, ranking_morris |
| `simulation_results_<TAG>.csv`                      | **All valid runs** with metrics + parameter values                       |
| `simulation_results_satisfactories_<TAG>.csv`       | **Only runs that passed the thresholds** (same columns as above)         |
| `inherited_ranges_<TAG>.csv`                        | Narrowed parameter ranges for the next iteration (auto-discovered)           |

The two `simulation_results_*` CSVs share the same schema and are both written when applicable. The choice of which one feeds `inherited_ranges_<TAG>.csv` is controlled by `cal.inheritance_source` (see Section 4 / Section 5).

#### `run_info.txt` — what's inside

Plain-text metadata summary written into every output folder. Designed for offline reproducibility: open it months later to recover exactly what was simulated without loading any `.rds`. Sections:

- **Header** — generation timestamp, output directory, R version, iteration label, `tag_prefix`
- **Simulation settings** — `project_path`, start/end dates, warmup years, threads, seed, outlet channel, temporal scale
- **Screening-specific** (Morris screening runs only): `m_traj`, `delta`, screening metrics, water-balance targets
- **Calibration-specific** (CAL runs only): `cal type` (soft/hard), `n_simulations`, `inheritance_source`, NSE/KGE/PBIAS thresholds, WB targets, `inherit from`, optional `inherit.drop` and `inherit.override`
- **Parameters used** — table of `short`, `par_name`, `min`, `max` for each parameter
- **Result summary** — `tag`, `n_simulations`, `n_valid`, `n_satisfactory`, `n_used_inheritance`, top-ranked parameter (screening), screening source folder (calibration)

Example excerpt:

```
==============================================================================
  swatplus-hmscal — CAL run summary
==============================================================================
  generated  : 2026-05-02 22:14:08 BRT
  output dir : C:/path/to/swatplus/txtinout/ITER3_v1
  config     : iteration3
  iteration  : 3
  tag_prefix : ITER3
------------------------------------------------------------------------------
Simulation settings
------------------------------------------------------------------------------
  start_date       : 1976-01-01
  end_date         : 1989-12-31
  warmup (years)   : 4
  threads          : 16
  seed             : 123
------------------------------------------------------------------------------
Calibration settings
------------------------------------------------------------------------------
  cal type            : hard
  n_simulations       : 250
  inheritance_source  : satisfactory
  threshold_nse        : 0.50
  threshold_kge        : 0.50
  threshold_abs_pbias  : 15
  WB diagnostic       : ET/P=0.48147, Q/P=0.51853, threshold_pct=15
  inherit from        : ITER2
------------------------------------------------------------------------------
Parameters used (18)
------------------------------------------------------------------------------
  short        par_name                                                  min          max
  cn2          cn2::cn2.hru | change = pctchg                       -25.0000      25.0000
  perco        perco::perco.hru | change = absval                     0.0000       1.0000
  ...
------------------------------------------------------------------------------
Result summary
------------------------------------------------------------------------------
  tag                    : ITER3_v1
  cal_type               : hard
  inheritance_source     : satisfactory
  n_simulations          : 250
  n_valid                : 248
  n_satisfactory         : 162
  n_used_inheritance     : 162
  morris_source             : C:/.../MORRIS3_v1
```

#### Calibration plots

| File                                            | Description                                                          |
|-------------------------------------------------|----------------------------------------------------------------------|
| `balance_scatter_<TAG>.tif`                     | ET/P vs Q_total/P scatter (ITER1 only)             |
| `performance_scatter_<TAG>.tif`                 | NSE vs KGE scatter (ITER2, ITER3, ITER4)         |
| `hydrograph_<TAG>.tif`                          | Observed vs simulated streamflow                                     |
| `params_boxplot_<TAG>.tif`                      | Normalized parameter distributions (accepted runs)                |

### Key Columns in Results CSV

**Iteration 1 (`type: "soft"`):**
- `rto_et`, `rto_wyld` — simulated ET/P and WYLD/P ratios
- `err_rto_et`, `err_rto_wyld` — relative errors vs targets
- `satisfactory_run` — TRUE if both errors <= threshold

**Iterations 2 to 4 (`type: "hard"`):**
- `NSE_mon` or `NSE_day` — Nash-Sutcliffe Efficiency
- `KGE_mon` or `KGE_day` — Kling-Gupta Efficiency
- `PBIAS_mon` or `PBIAS_day` — Percent Bias
- `obj_total` — Combined objective function
- `satisfactory_run` — TRUE if NSE, KGE, and PBIAS pass thresholds

### Loading Saved Results in R

```r
# Load full results object (filename uses the TAG of the producing folder)
res <- readRDS("C:/path/to/your/swatplus/txtinout/ITER2_v1/simulation_results_ITER2_v1.rds")

# Access components
res$tag                  # "ITER2_v1" — self-identifying
res$config               # YAML config used
res$param_info           # Parameters and ranges
res$ranking_morris       # Morris ranking that fed the narrowing
res$results              # All valid runs with metrics
res$satisfactory_results # accepted runs only
res$sim                  # Raw SWATrunR simulation output

# Quick summary
cat("Total runs:", nrow(res$results), "\n")
cat("Satisfactory:", nrow(res$satisfactory_results), "\n")
```

For Morris screening results the equivalent file is `morris_results_<TAG>.rds`.

---

## 10. Adapting to Your Catchment

### Step-by-Step Customization

#### 1. Set your project path

Edit all three YAML configs:
```yaml
project_path: "C:/path/to/your/swatplus/txtinout"
```

#### 2. Set your outlet channel unit

In the YAML config files, find the output definitions and update the `unit` to match your outlet channel number:
```yaml
morris_outputs:
  - file: "channel_sd"
    variable: "flo_out"
    unit: 33                  # <-- Change to YOUR outlet channel number
```

> **How to find your outlet channel:** Open `channel-lte.cha` in your TxtInOut folder.
> The outlet is typically the channel with the largest drainage area or the one
> corresponding to your gauging station location.

#### 3. Update observed data file names

```yaml
obs_files:
  monthly: "obs_flow_monthly.csv"    # Your monthly observed streamflow file
  daily:   "obs_flow_daily.csv"      # Your daily observed streamflow file
```

#### 4. Set simulation period and warmup

```yaml
start_date: "1976-01-01"     # Simulation start (including warmup)
end_date:   "1989-12-31"     # Simulation end
warmup: 4                    # Years to skip (warmup period)
obs_start_date: "1980-01-01" # First date of observed data to use
```

#### 5. Set water balance targets (Iteration 1)

Calculate from observed data or literature for your catchment:
```yaml
cal:
  targets:
    et_rto: 0.48147          # Long-term ET / Precipitation ratio
    wyld_rto: 0.51853        # Long-term Water Yield / Precipitation ratio (= 1 - ET/P)
```

> **Tip:** `ET/P + WYLD/P ≈ 1.0`. Use long-term annual averages from your catchment.
> ET can be estimated from remote sensing (MODIS, GLEAM) or water balance closure.

#### 6. Adjust parameters for your region

Iteration 1 parameters are defined directly in `config/iteration1.yaml`. You can:
- **Add parameters:** Add new entries to the `parameters` list
- **Remove parameters:** Delete entries you don't need
- **Change ranges:** Adjust `min` and `max` values

For Iterations 2 to 4, the `param_map` section controls which parameters are carried forward.

#### 7. Adjust calibration thresholds

Based on your data quality and expectations:
```yaml
# Iteration 2 (monthly) — stricter thresholds
cal:
  threshold_nse: 0.70
  threshold_kge: 0.70
  threshold_abs_pbias: 10

# Iteration 3 (daily) — relaxed thresholds
cal:
  threshold_nse: 0.50
  threshold_kge: 0.50
  threshold_abs_pbias: 15
```

> **Reference:** Thresholds follow Moriasi et al. (2015) guidelines.
> NSE > 0.50 is "satisfactory", > 0.65 is "good", > 0.75 is "very good".

#### 8. Adjust computational resources

```yaml
threads: 16                  # Number of CPU cores for parallel SWAT+ runs (range: check your devide core number)
cal:
  n_simulations: 250         # More simulations = better exploration but slower (range: unlimited)
morris:
  m_traj: 15                 # More trajectories = more robust Morris indices (range: 10-20 traj)
```

---

## 11. Troubleshooting

### Common Errors

**"Error: Usage: Rscript run_morris.R config/iterationN.yaml"**
- You forgot to set `ITERATION_CONFIG` before `source()`, or didn't pass the YAML path as argument.

**"Error in yaml::read_yaml(...) : file not found"**
- Check your working directory: `getwd()`. It must be the `swatplus-hmscal/` folder.

**"Error: ABORTING: some parameters exceed the cal_parms.cal bounds"**
- Iteration 1 range check found parameters that would produce values outside SWAT+ limits.
- Review the printed table and adjust `min`/`max` in `config/iteration1.yaml`.

**"Error: No folder found with prefix '...'"**
- The inheritance chain is broken. For Iteration 2, you need Iteration 1 CAL results first.
- Run the previous iteration before proceeding.

**"Error: Column 'Flow' not found" or empty observed data**
- Check your observed CSV file: must have columns named exactly `Date` and `Flow`.
- Check date format: must be `YYYY-MM-DD`.
- Check `obs_start_date` falls within your observed data range.

**"Satisfactory runs: 0"**
- No runs passed the thresholds. Options:
  - Relax thresholds (lower NSE/KGE, increase PBIAS tolerance)
  - Widen parameter ranges
  - Increase `n_simulations` for better coverage
  - Check if observed data is correctly aligned with simulated period

**"Error in run_swatplus(...)" or SWAT+ crashes**
- Test your SWAT+ model independently first (see Section 3).
- Check that `project_path` points to a valid TxtInOut folder.
- Reduce `n_thread` if you're running out of memory.

**"Can't subset columns past the end. Locations 50 and 51 don't exist." (lsunit_wb_aa.txt)**
- Triggered when an output reads `lsunit_wb_aa.txt` on SWAT+ rev. ≥ 2026.61.0.2.61.
- Root cause: the SWAT+ writer emits 51 names in the header but only 49 numeric
  columns per data row (`plant_cov` and `mgt_ops` are in the header but never
  written). SWATrunR 0.9.4 `read_output_i()` does `data[, 1:length(col_names_i)]`
  and crashes — losing the whole `sim` object even though the SWAT+ runs
  themselves succeeded.
- **Fix: rebuild SWATrunR from source with a 6-line guard in `read_output_i()`.**
  An in-memory monkey-patch via `utils::assignInNamespace()` does **not** work
  here, because `run_swatplus()` calls `read_output_i()` inside `foreach %dopar%`
  workers — separate R processes that load SWATrunR fresh from disk and don't
  see master-side namespace edits. The fix has to be on disk.
  Step-by-step instructions, the exact source patch, and propagation to a
  remote machine are documented in the sister project's README:
  [`../swatplus-calibrator/README.md`](../swatplus-calibrator/README.md) →
  section "Known issues and workarounds". One-time install per machine; reverts
  with `remotes::install_github("chrisschuerz/SWATrunR", upgrade = "never")`.
- Validated against SWATrunR 0.9.4 (SHA `82e6353001480027303b230adf57bf5c7ba7057e`),
  R 4.3.3, SWAT+ rev. 2026.61.0.2.61.

### Performance Tips

- **The Iteration 1 screening** is fast (few runs, annual output). Start here to verify the setup works.
- **Iteration 2-3 daily simulations** can be slow. Start with fewer simulations (e.g., 100) to test, then increase.
- **Disk space:** SWATrunR creates temporary run folders. The scripts clean them before each run, but check disk space if runs fail mid-execution.
- **Memory:** Large numbers of simulations can consume significant RAM. Monitor with Task Manager.

### Re-running a Step

Simply run the same step again — a new versioned folder (v2, v3...) will be created automatically. Previous results are preserved.

```r
# Re-run MORRIS2 (will create MORRIS2_v2/, then v3, ...)
ITERATION_CONFIG <- "config/iteration2.yaml"
source("run_morris.R")
```

> **Note:** When Iteration 2+ inherits ranges, `find_inherited_ranges()` always picks the **latest version** folder matching the source TAG. If you re-run ITER1, MORRIS2 will automatically pick up the newest `ITER1_v<N>/inherited_ranges_ITER1_v<N>.csv`.

---

## Quick Reference

```r
# ---- Full pipeline (7 steps) ----
source("master.R")

# ---- Resume from step 3 ----
START_FROM <- 3
source("master.R")

# ---- Run a single Morris screening ----
ITERATION_CONFIG <- "config/iteration1.yaml"
source("run_morris.R")

# ---- Run single calibration ----
ITERATION_CONFIG <- "config/iteration2.yaml"
source("run_cal.R")

# ---- Run ITER4 (depth refinement of ITER3) ----
ITERATION_CONFIG <- "config/iteration4.yaml"
source("run_cal.R")

# ---- From command line ----
# Rscript run_morris.R  config/iteration1.yaml
# Rscript run_cal.R  config/iteration1.yaml
# Rscript master.R

# ---- Load saved results ----
res <- readRDS("C:/path/to/your/swatplus/txtinout/ITER2_v1/simulation_results_ITER2_v1.rds")
```

---

## Release history

| Version | Date       | Notes                                                                                             |
|---------|------------|---------------------------------------------------------------------------------------------------|
| 1.1.0   | 2026-10-27 | Renamed to **HMSCal** (`swatplus-hmscal`). Adds the split-sample validation workflow (`run_val.R`, `config/validation.yaml`), the Eckhardt BFI module (`R/bfi_eckhardt.R`), baseflow diagnostics, the Iteration 3 Iteration 2 config, and worked example observed data. This is the version used for the manuscript results. |
| 1.0.0   | 2026-04-26 | First archived release, as `swatplus-hms-calibration`. Calibration cascade only, without the validation workflow or the BFI modules. |

---

## Citation

If you use this software, please cite it via the metadata in [`CITATION.cff`](CITATION.cff). After a tagged release, a versioned DOI is automatically minted by Zenodo and shown as a badge at the top of this README.

GitHub natively renders `CITATION.cff` as a "Cite this repository" button on the repository landing page (top right), with ready-to-paste BibTeX/APA strings.

The badge above carries the **concept DOI**, which always resolves to the latest version. To cite the exact code behind a given result, use the versioned DOI of the corresponding release rather than the concept DOI. Note that the v1.0.0 deposit was archived under the former name, `swatplus-hms-calibration`.

---

## Contact

For technical questions, bug reports, or feature requests, please open an issue on GitHub: <https://github.com/ECRIPPEL/swatplus-hmscal/issues>.

Author: Elzon Cassio Rippel — ORCID [0000-0002-8391-4435](https://orcid.org/0000-0002-8391-4435).

---

## License

Released under the MIT License — see [`LICENSE`](LICENSE) for the full text.

This tool depends on [SWAT+](https://swat.tamu.edu/software/plus/) and [SWATrunR](https://chrisschuerz.github.io/SWATrunR/), which are distributed under their own licenses.
