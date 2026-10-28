# SWATrunR_mod for SWAT+ rev. 2026.61.0.2.61

HMSCal runs on top of [SWATrunR](https://github.com/chrisschuerz/SWATrunR) (Schürz, GPL-3), which handles the SWAT+ executable calls, parameter writing and parallel execution. A stock install of SWATrunR 0.9.4 will **not** get you through the HMSCal pipeline on a current R toolchain. Two separate problems get in the way, and this folder documents both along with the fixes.

Read this before your first run. The second problem in particular does not fail fast: it wastes the entire simulation batch and only then aborts.

| Problem | Origin | Where it bites | How it is fixed here |
|---|---|---|---|
| `readr::read_table2()` no longer exists | R ecosystem | Package installation | Pin `readr` below 2.2.0 |
| `lsunit_wb_aa.txt` header is wider than its data rows | SWAT+ writer | After every SWAT+ thread finishes | Source patch, rebuild from source |

### Why a frozen reader meets a moving model

Neither problem is really a bug in SWATrunR's logic. Both are drift. SWATrunR 0.9.4 is the last upstream state (commit `82e6353`), and it has stood still while everything around it moved: `readr` dropped a function it calls, and SWAT+ kept shipping revisions whose output files it has to parse. A reader that is pinned while the writer is updated continuously will eventually read a file whose shape it was never told about. That is exactly what `lsunit_wb_aa.txt` is.

This is why the SWAT+ revision is in the title of this document. The column mismatch is a property of what **that** revision's Fortran writer emits, not a universal truth about SWAT+.

| Component | Version used here | Notes |
|---|---|---|
| SWAT+ | rev. 2026.61.0.2.61 | The revision whose `lsunit_wb_aa.txt` writer motivated the guard |
| SWATrunR | 0.9.4, commit `82e6353` | Upstream HEAD at the time; the patch is a diff against this exact commit |
| R | 4.3.3 | |
| readr | 2.1.5 | Any 2.1.x works. 2.2.0 and later do not |

**If you are on a different SWAT+ revision.** The readr problem is independent of SWAT+ and will affect you regardless. The guard is written defensively rather than against specific numbers: it hardcodes neither `51`, nor `49`, nor the file name. If a later revision fixes the writer, the guard simply never fires and costs you nothing. If a later revision introduces a mismatch in some other output file, the guard catches that one too and names the columns it dropped. So the patch is useful beyond the revision in the title, but only the revision in the title was actually tested.

---

## Install it

From R, with the repository root as your working directory:

```r
source("SWATrunR_mod/install_swatrunr_patched.R")
install_swatrunr_patched()
```

That is the whole thing. The script downloads SWATrunR 0.9.4 at the pinned upstream commit, pins `readr` if yours is too new, applies the guard, installs, and then verifies the **installed** package actually carries the fix. Restart R afterwards.

You need R (tested on 4.3.3) and an internet connection. You do **not** need git, RTools, or any compiler: SWATrunR has no `src/` directory, so nothing is compiled.

Expected output:

```
[swatrunr-patch] readr 2.1.5 is below 2.2.0, keeping it
[swatrunr-patch] downloading SWATrunR at 82e6353
[swatrunr-patch] guard inserted at line 104 of R/swatplus_handle_output.R
[swatrunr-patch] installing into C:/Users/you/AppData/Local/R/win-library/4.3
...
[swatrunr-patch] OK. SWATrunR 0.9.4 installed with the header guard.
```

If it stops with an error it will say which step failed. The script refuses to report success unless the guard is present in the installed package, so a silent half-install is not possible. To install somewhere other than your first library, pass `install_swatrunr_patched(lib = "path/to/lib")`.

The rest of this document explains what the two problems are, why the fix has to be compiled into the package, and how to do the same thing by hand.

---

## Problem 1: installation fails on readr 2.2.0 and newer

### Symptom

Installing SWATrunR 0.9.4 with `remotes::install_github()` on a clean R 4.3.3 fails:

```
Erro: objeto 'read_table2' não foi exportado por 'namespace:readr'
Execução interrompida
ERROR: lazy loading failed for package 'SWATrunR'
```

### Cause

SWATrunR calls `readr::read_table2()` in `R/swatplus_parameter_write.R`. That function was soft-deprecated in readr 2.0.0 (July 2021, when readr was rewritten on top of `vroom`) and removed outright in readr 2.2.0. Since current CRAN readr is at or above 2.2.0, a default install breaks before the package ever loads.

### Fix

Pin readr below 2.2.0 before installing SWATrunR:

```r
remotes::install_version("readr", version = "2.1.5", upgrade = "never")
remotes::install_github("chrisschuerz/SWATrunR", upgrade = "never")
```

Any readr in the 2.1.x line works, because `read_table2()` is still exported there (deprecated, but present). The HMSCal results were produced with readr 2.1.5.

Keep readr pinned. A later `update.packages()` will pull readr past 2.2.0 and break SWATrunR again. If you prefer not to pin, replacing `read_table2` with `read_table` throughout the SWATrunR source is behaviour-preserving, since `read_table` in readr 2.0 and later has the semantics of the old `read_table2`. We did not take that route, so it is untested here.

---

## Problem 2: `lsunit_wb_aa.txt` aborts the run after it has already finished

This is the expensive one.

### Symptom

Any call to `run_swatplus()` that requests an output from `lsunit_wb_aa` fails during post-processing:

```
Error in `.[, 1:length(col_names_i)]`:
! Can't subset columns past the end.
ℹ Locations 50 and 51 don't exist.
ℹ There are only 49 columns.
```

The crucial detail is the timing. SWAT+ itself runs correctly, every parallel thread completes, and the failure happens while SWATrunR is assembling the results. The `sim` object is destroyed on its way out. On a multi-member ensemble over a long simulation period that is tens of minutes to hours of CPU time lost with nothing written to disk.

Reproduced with SWAT+ rev. 2026.61.0.2.61 and triggered by any `define_output(file = "lsunit_wb_aa", ...)`.

### Cause

Two independent defects meet.

**The SWAT+ writer emits a malformed file.** `lsunit_wb_aa.txt` (annual mean water balance per landscape unit) is written with 51 names in the header line but only 49 numeric columns per data row. `plant_cov` and `mgt_ops` appear in the header and are never written to the rows.

**SWATrunR trusts the header.** `read_output_i()` in `R/swatplus_handle_output.R` does:

```r
fread(thread_path %//% output_i$file_full[1], skip = 3) %>%
  as_tibble(.) %>%
  .[, 1:length(col_names_i)] %>%
  ...
```

When `col_names_i` has 51 entries and the data has 49 columns, tibble's strict subsetting throws, the error propagates up, and `run_swatplus()` returns nothing.

### Why a runtime patch does not work

The obvious workaround is to override the function after loading the package:

```r
utils::assignInNamespace("read_output_i", patched, ns = "SWATrunR")   # does NOT work
```

This passes synthetic tests on the master process and then fails in production. `run_swatplus()` calls `read_output_i()` inside `foreach %dopar%` workers spawned by `parallel::makeCluster()`. PSOCK workers are fresh R processes, and when functions are serialized to them namespaces travel **by name**, not by value. Each worker calls `loadNamespace("SWATrunR")` and reads the package from disk. Edits made in the master's memory never reach them.

Confirmed empirically on an 86-simulation run: the patch banner printed in the master log, all 86 simulations completed, and a worker then raised the identical "Can't subset columns past the end" error.

The practical rule: anything that has to survive `foreach %dopar%` in SWATrunR must be on disk, in the package binary itself.

### The fix

A guard in `read_output_i()` that truncates the header to the width the data actually has, and warns about what it dropped. Shipped here as [`01-read_output_i-header-guard.patch`](01-read_output_i-header-guard.patch), a diff against `R/swatplus_handle_output.R` at upstream commit `82e6353` (v0.9.4):

```r
n_data <- ncol(raw_data)
if (length(col_names_i) > n_data) {
  dropped <- col_names_i[(n_data + 1):length(col_names_i)]
  warning(sprintf(
    "[SWATrunR] '%s': header has %d cols but data has %d. Truncating header (dropping: %s).",
    output_i$file_full[1], length(col_names_i), n_data,
    paste(dropped, collapse = ", ")), call. = FALSE)
  col_names_i <- col_names_i[seq_len(n_data)]
}
```

It fixes `lsunit_wb_aa.txt` but will equally absorb any other SWAT+ output whose header is wider than its rows, and it tells you which names it discarded instead of failing silently.

The patch was validated against a synthetic file with a 51-name header and 49 data columns, in the master process and inside doSNOW workers, which is the configuration that actually failed before.

---

## Installing by hand

Use this if you want to inspect every step, or if the script cannot reach GitHub from your network. Needs git in addition to R. Still no compiler.

```bash
# 1. Pin readr first (see Problem 1), from R:
#    remotes::install_version("readr", version = "2.1.5", upgrade = "never")

# 2. Clone upstream at the exact commit the patch applies to
git clone https://github.com/chrisschuerz/SWATrunR.git /c/temp/SWATrunR-src
cd /c/temp/SWATrunR-src
git checkout 82e6353001480027303b230adf57bf5c7ba7057e

# 3. Apply the patch shipped in this folder
git apply /path/to/swatplus-hmscal/SWATrunR_mod/01-read_output_i-header-guard.patch

# 4. Build and install over any existing SWATrunR
"/c/Program Files/R/R-4.3.3/bin/R.exe" CMD INSTALL /c/temp/SWATrunR-src/
```

Step 4 can equally be done from inside R with `install.packages("C:/temp/SWATrunR-src", repos = NULL, type = "source")`.

Use `git apply` from inside the clone, as above. The patch is a git diff with LF line endings, and on Windows git checks the file out as CRLF, so the GNU `patch` utility rejects it with "different line endings". `git apply` resolves this against the index and applies cleanly. If you must use `patch`, convert first with `dos2unix` or pass `--binary`. The installer script sidesteps this entirely by editing the source lines directly rather than applying a diff.

### Verifying

From R, after restarting the session:

```r
any(grepl("Truncating header", deparse(SWATrunR:::read_output_i)))
# must return TRUE
```

If it returns `FALSE` you are still on a stock install and Problem 2 will hit you.

### Copying to a second machine

If both machines run the same R version, you do not need RTools on the second one. Run `find.package("SWATrunR")` on the patched machine, zip that whole folder, and unzip it into the `.libPaths()` location of the target (rename the existing `SWATrunR/` folder first, as a backup). Restart R and run the verification above.

### Reverting to stock

```r
remotes::install_github("chrisschuerz/SWATrunR", upgrade = "never")
```

---

## Scope and attribution

SWATrunR is the work of Christoph Schürz and is distributed under **GPL-3**. HMSCal is MIT. For that reason this folder ships only the patch, not a copy of the package: you fetch SWATrunR from its own repository, under its own licence, and apply a diff. Nothing here redistributes upstream code.

The patch touches one function in one file. Everything else in SWATrunR is untouched, and the SWAT+ executable interface is not involved, so model output is unaffected. The guard changes only how a malformed output file is parsed.

The underlying defect in Problem 2 is in the SWAT+ Fortran writer, not in SWATrunR. The R-side guard keeps the run alive; it does not recover the two columns that SWAT+ never wrote. If you need `plant_cov` or `mgt_ops` from `lsunit_wb_aa`, neither this patch nor stock SWATrunR can give them to you.
