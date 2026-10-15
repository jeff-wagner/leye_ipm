# =============================================================================
# Recruitment prior sensitivity runs
#
# Refits the IPM with ONLY the recruitment hyperparameters changed. Everything
# else -- data, structure, sampler settings -- is whatever scripts/LEYE_IPM.R
# currently defines, so these fits stay comparable to the main one by
# construction rather than by anyone remembering to keep them in step.
#
# SETTINGS
#   cohortmean  pi centred on the MEAN ACROSS modern cohorts rather than on the
#               typical cohort. LEYE_recruitment_two_era.R fits
#                   logit(pi_c) = mu_pi + b_era*modern + eps_c
#               and reports pi_M = ilogit(mu_pi + b_era), i.e. the cohort with
#               its random effect at zero. But the IPM uses pi in an EXPECTED
#               count, R = N*(f/2)*psi, with no cohort random effect of its own,
#               so arguably it wants E[pi_c]. With sd_pi ~ 1.39 the two differ by
#               about a third (0.054 vs 0.072) through Jensen's inequality.
#               Hold this one loosely: sd_pi rests on 13 cohorts and 18 returning
#               birds, and the 1997 cohort alone sits at pi_c = 0.58 against
#               0.03-0.18 elsewhere. Since sd_pi is shared across eras, that one
#               cohort inflates the gap for the modern cohorts too.
#   vague       Beta(1,1) on both recruitment parameters: what the counts alone
#               would say. This used to be produced by hand-editing constants,
#               which is why it had no script and no record of how it was made.
#
# USAGE
#   Rscript scripts/LEYE_IPM_sensitivity.R              # both settings
#   Rscript scripts/LEYE_IPM_sensitivity.R cohortmean   # just one
#
# Outputs are deliberately NOT the canonical cache name: the report and the
# tLTRE both read LEYE_IPM_samples.rds, and a companion fit must never be
# mistaken for the main one.
# =============================================================================

library(nimble)
library(parallel)
library(coda)

SKIP_RUN <- TRUE
source("scripts/LEYE_IPM.R") # ipmCode, constants, data, inits, monitors

## ---------------------------------------------------------------------------
## Settings
## ---------------------------------------------------------------------------
# Beta fit to the cohort-averaged modern pi, by method of moments.
cohortmean_prior <- function() {
  te <- do.call(rbind, readRDS("LEYE_recruitment_two_era.rds"))
  pic <- grep("^pi_c\\[", colnames(te), value = TRUE)
  N_COH_90 <- 5 # 1995-1999; modern cohorts follow (see `modern` in that module)
  modern <- pic[(N_COH_90 + 1):length(pic)]
  stopifnot(length(modern) > 0)
  x <- rowMeans(te[, modern, drop = FALSE])
  mu <- mean(x)
  k <- mu * (1 - mu) / var(x) - 1
  list(a_pi = mu * k, b_pi = (1 - mu) * k)
}

SETTINGS <- list(
  cohortmean = list(
    out = "LEYE_IPM_samples_cohortmean.rds",
    label = "cohort-averaged pi",
    make = cohortmean_prior # kappa unchanged: it is a scalar, no Jensen gap
  ),
  vague = list(
    out = "LEYE_IPM_samples_vaguepriors.rds",
    label = "vague Beta(1,1) on both",
    make = function() list(a_pi = 1, b_pi = 1, a_kap = 1, b_kap = 1)
  )
)

args <- commandArgs(trailingOnly = TRUE)
which_run <- if (length(args)) args else names(SETTINGS)
unknown <- setdiff(which_run, names(SETTINGS))
if (length(unknown)) {
  stop("Unknown setting(s): ", paste(unknown, collapse = ", "),
       ". Available: ", paste(names(SETTINGS), collapse = ", "))
}

## ---------------------------------------------------------------------------
## Runner -- same sampler settings as scripts/LEYE_IPM_run_parallel.R
## ---------------------------------------------------------------------------
nchains <- 3
seeds <- 1:nchains
niter <- 100000
nburnin <- 30000
thin <- 10

run_one <- function(key) {
  s <- SETTINGS[[key]]
  cs <- constants
  for (nm in names(s$make())) cs[[nm]] <- s$make()[[nm]]

  cat(sprintf(
    "\n=== %s (%s) ===\n  pi_rec ~ Beta(%.4f, %.4f)   kappa ~ Beta(%.4f, %.4f)\n",
    key, s$label, cs$a_pi, cs$b_pi, cs$a_kap, cs$b_kap
  ))
  cat(sprintf("  main fit uses Beta(%.4f, %.4f) and Beta(%.4f, %.4f)\n",
              constants$a_pi, constants$b_pi, constants$a_kap, constants$b_kap))

  inits_per_chain <- lapply(seeds, function(x) {
    set.seed(x)
    inits()
  })

  # Stop the cluster exactly once, on exit. Calling stopCluster() explicitly as
  # well as registering it here closes the sockets twice: the first setting
  # finishes and saves, then the next makeCluster() inherits broken state and
  # parLapply dies with "invalid connection" partway through. on.exit also
  # covers the error path, which an explicit call at the end does not.
  cl <- makeCluster(min(nchains, detectCores() - 1))
  on.exit(stopCluster(cl), add = TRUE)
  t0 <- Sys.time()
  chains <- parLapply(
    cl, seq_len(nchains),
    function(i, seeds, ipmCode, constants, data, inits_per_chain, monitors,
             niter, nburnin, thin) {
      library(nimble)
      set.seed(seeds[i])
      Rmodel <- nimbleModel(ipmCode,
        constants = constants, data = data,
        inits = inits_per_chain[[i]], calculate = FALSE
      )
      Rmcmc <- buildMCMC(configureMCMC(Rmodel, monitors = monitors))
      Cmodel <- compileNimble(Rmodel)
      Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
      runMCMC(Cmcmc,
        niter = niter, nburnin = nburnin, thin = thin,
        nchains = 1, setSeed = seeds[i], samplesAsCodaMCMC = TRUE
      )
    },
    seeds = seeds, ipmCode = ipmCode, constants = cs, data = data,
    inits_per_chain = inits_per_chain, monitors = monitors,
    niter = niter, nburnin = nburnin, thin = thin
  )

  samples <- as.mcmc.list(chains) |>
    lapply(function(x) replace(x, is.na(x), 0)) |>
    as.mcmc.list()
  saveRDS(stamp_ipm_cache(samples), s$out)
  cat(sprintf(
    "  done in %.1f min -> %s\n",
    as.numeric(difftime(Sys.time(), t0, units = "mins")), s$out
  ))
}

for (k in which_run) run_one(k)
cat("\nAll requested sensitivity fits complete.\n")
