# =============================================================================
# Parallel MCMC for the LEYE IPM -- one chain per core.
#
# WHY a worker function: NIMBLE's compiled model/MCMC are bound to the R process
# that built them and CANNOT be exported to workers. So each worker must do the
# full build itself: nimbleModel -> configureMCMC -> buildMCMC -> compileNimble
# -> runMCMC. We pass only plain R data (constants/data/inits) into the workers.
#
# Output: a coda mcmc.list identical in form to runMCMC(..., nchains=k), so all
# your existing summary/plot scripts work unchanged on `samples`.
# =============================================================================

library(parallel)

## ---------------------------------------------------------------------------
## 1. Build the data objects ONCE in the main session (cheap, no NIMBLE compile)
##    Source ONLY the data-prep portion of your IPM script. Easiest: refactor
##    Sections 1-3 of leye_ipm_full_covariates.R into a function that returns
##    the model pieces. Below assumes you've sourced that script up through the
##    definition of: ipmCode, constants, data, inits (function), monitors.
## ---------------------------------------------------------------------------
# Option A (simplest): source the prep but DO NOT run Section 4/5.
#   Wrap Sections 4-5 of the IPM script in `if (!exists('SKIP_RUN'))` and set
#   SKIP_RUN <- TRUE before sourcing, so sourcing only defines the objects.
SKIP_RUN <- TRUE
source("scripts/LEYE_IPM.R") # defines ipmCode, constants, data, inits, monitors
# (ipm_years, T, cjs, etc. also become available for post-processing)

## ---------------------------------------------------------------------------
## 2. Worker function: builds + runs ONE chain. Everything NIMBLE happens here.
## ---------------------------------------------------------------------------
run_one_chain <- function(
  seed,
  ipmCode,
  constants,
  data,
  inits_list,
  monitors,
  niter,
  nburnin,
  thin
) {
  library(nimble)
  set.seed(seed)

  Rmodel <- nimbleModel(
    ipmCode,
    constants = constants,
    data = data,
    inits = inits_list,
    calculate = FALSE
  )
  conf <- configureMCMC(Rmodel, monitors = monitors)
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

  out <- runMCMC(
    Cmcmc,
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    nchains = 1,
    setSeed = seed,
    samplesAsCodaMCMC = TRUE
  )
  return(out) # a single coda mcmc object
}

## ---------------------------------------------------------------------------
## 3. Launch: one chain per core
## ---------------------------------------------------------------------------
nchains <- 3
seeds <- 1:nchains
niter <- 100000
nburnin <- 30000
thin <- 10

# Generate DISTINCT inits per chain (inits is a function in the IPM script).
# Distinct starting values across chains are what make Rhat meaningful.
inits_per_chain <- lapply(seeds, function(s) {
  set.seed(s)
  inits()
})

ncores <- min(nchains, detectCores() - 1)
cl <- makeCluster(ncores)
on.exit(stopCluster(cl), add = TRUE) # ensure cleanup even on error

t0 <- Sys.time()
# Pass EVERYTHING the worker needs as explicit arguments. parLapply forwards the
# trailing args to the function for every element, so each worker gets its own
# copies of ipmCode/constants/data/monitors/settings. Per-chain values (seed,
# inits) are selected by the index i inside the worker. No clusterExport needed.
chain_list <- parLapply(
  cl,
  seq_len(nchains),
  function(
    i,
    seeds,
    ipmCode,
    constants,
    data,
    inits_per_chain,
    monitors,
    niter,
    nburnin,
    thin
  ) {
    library(nimble)
    set.seed(seeds[i])
    Rmodel <- nimbleModel(
      ipmCode,
      constants = constants,
      data = data,
      inits = inits_per_chain[[i]],
      calculate = FALSE
    )
    conf <- configureMCMC(Rmodel, monitors = monitors)
    Rmcmc <- buildMCMC(conf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    runMCMC(
      Cmcmc,
      niter = niter,
      nburnin = nburnin,
      thin = thin,
      nchains = 1,
      setSeed = seeds[i],
      samplesAsCodaMCMC = TRUE
    )
  },
  seeds = seeds,
  ipmCode = ipmCode,
  constants = constants,
  data = data,
  inits_per_chain = inits_per_chain,
  monitors = monitors,
  niter = niter,
  nburnin = nburnin,
  thin = thin
)
stopCluster(cl)
cat(sprintf(
  "Parallel run done in %.1f min (%d chains on %d cores).\n",
  as.numeric(difftime(Sys.time(), t0, units = "mins")),
  nchains,
  ncores
))

## ---------------------------------------------------------------------------
## 4. Recombine into a coda mcmc.list -> drop-in replacement for `samples`
## ---------------------------------------------------------------------------
library(coda)
library(tidyr)
samples <- as.mcmc.list(chain_list) |>
  lapply(function(x) replace(x, is.na(x), 0))

# Quick convergence peek (your full summary script can run as-is after this)
cat("\nTop Rhat values:\n")
gd <- gelman.diag(samples, multivariate = FALSE)$psrf
print(round(head(gd[order(-gd[, 1]), ], 10), 3))
cat(sprintf(
  "\nParams with Rhat > 1.05: %d\n",
  sum(gd[, 1] > 1.05, na.rm = TRUE)
))

saveRDS(samples, "LEYE_IPM_samples_20260803.rds") # cache to avoid refitting
