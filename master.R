# ==============================================================================
# master.R — Orchestrator for the four HMSCal iterations.
#
# A Morris screening precedes Iterations 1 to 3; Iteration 4 reuses the
# classification from the third screening. That is 7 steps in total:
#
#   MORRIS1 -> ITER1 -> MORRIS2 -> ITER2 -> MORRIS3 -> ITER3 -> ITER4
#
# Runs all 7 steps sequentially.
# Can be interrupted and resumed from any step.
#
# USAGE:
#   source("master.R")                    # run all steps
#   START_FROM <- 3; source("master.R")   # resume from step 3 (MORRIS2)
# ==============================================================================

# Starting step (default = 1)
if (!exists("START_FROM")) START_FROM <- 1

steps <- list(
  list(step = 1, label = "MORRIS1 - Morris screening, annual scale",
                                                 script = "run_morris.R", config = "config/iteration1.yaml"),
  list(step = 2, label = "ITER1 - Iteration 1, annual water balance (RC + BFI proxy)",
                                                 script = "run_cal.R", config = "config/iteration1.yaml"),
  list(step = 3, label = "MORRIS2 - Morris screening, monthly scale",
                                                 script = "run_morris.R", config = "config/iteration2.yaml"),
  list(step = 4, label = "ITER2 - Iteration 2, monthly streamflow + RC + BFI proxy",
                                                 script = "run_cal.R", config = "config/iteration2.yaml"),
  list(step = 5, label = "MORRIS3 - Morris screening, daily scale",
                                                 script = "run_morris.R", config = "config/iteration3.yaml"),
  list(step = 6, label = "ITER3 - Iteration 3, daily streamflow + RC + BFI proxy (routing params enter)",
                                                 script = "run_cal.R", config = "config/iteration3.yaml"),
  list(step = 7, label = "ITER4 - Iteration 4, same daily criteria on the narrowed ranges",
                                                 script = "run_cal.R", config = "config/iteration4.yaml")
)

cat("==============================================================\n")
cat(" HMSCal - HIERARCHICAL MULTISCALE CALIBRATION OF SWAT+\n")
cat(" MORRIS1 -> ITER1 -> MORRIS2 -> ITER2 -> MORRIS3 -> ITER3 -> ITER4\n")
cat("==============================================================\n")
cat(" Starting from step:", START_FROM, "\n\n")

for (s in steps) {
  if (s$step < START_FROM) {
    cat("[SKIP] Step", s$step, "-", s$label, "\n")
    next
  }

  cat("\n==============================================================\n")
  cat(" STEP", s$step, "/", length(steps), "-", s$label, "\n")
  cat("==============================================================\n")

  t_start <- Sys.time()

  ITERATION_CONFIG <- s$config
  .MASTER_STATE <- list(steps = steps, s = s, START_FROM = START_FROM, t_start = t_start)
  source(s$script, local = FALSE)
  list2env(.MASTER_STATE, envir = globalenv())
  rm(.MASTER_STATE)

  t_end <- Sys.time()
  elapsed <- round(difftime(t_end, t_start, units = "mins"), 1)

  cat("\n[DONE] Step", s$step, "-", s$label, "in", elapsed, "minutes\n")
}

cat("\n==============================================================\n")
cat(" CALIBRATION COMPLETE!\n")
cat("==============================================================\n")
