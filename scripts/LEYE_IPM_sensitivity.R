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
#               about a third (0.054 vs 0.072): eps_c is symmetric about zero on
#               the LOGIT scale, but ilogit bends, so back-transforming the
#               typical cohort is not the same as averaging the back-transformed
#               cohorts (Jensen's inequality -- averaging before and after a
#               curved transformation gives different answers).
#               Hold this one loosely: sd_pi rests on 13 cohorts and 18 returning
#               birds, and the 1997 cohort alone sits at pi_c = 0.58 against
#               0.03-0.18 elsewhere. Since sd_pi is shared across eras, that one
#               cohort inflates the gap for the modern cohorts too.
#   drop97      both recruitment parameters re-derived from the two-era fit with
#               the 1997 cohort REMOVED (scripts/LEYE_recruitment_drop_cohort.R).
#               That cohort returned 7 of 11 chicks -- pi_c = 0.58 against <=0.12
#               for every other cohort in either era -- and it sets the 1990s
#               baseline the era contrast is measured against. Dropping it moves
#               pi from 0.054 to 0.044 and, unlike the other settings, also makes
#               the prior TIGHTER (sd 0.021 vs 0.038). The others all shift the
#               mean and leave the width alone, so this is the one direction the
#               sensitivity had not yet probed: a harder-pinned pi could shift
#               the recruitment/immigration split even with lambda unmoved.
#               kappa moves too (0.47 -> 0.40), so both come from the same fit
#               rather than mixing eras of evidence.
#   vague       Beta(1,1) on both recruitment parameters: what the counts alone
#               would say. This used to be produced by hand-editing constants,
#               which is why it had no script and no record of how it was made.
#
# USAGE
#   Rscript scripts/LEYE_IPM_sensitivity.R              # every setting
#   Rscript scripts/LEYE_IPM_sensitivity.R cohortmean   # just one
#   Rscript scripts/LEYE_IPM_sensitivity.R restamp      # fix provenance, no refit
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
# Method-of-moments Beta for a vector of posterior draws. This is the formula
# that produced the literals in LEYE_IPM.R -- moment-matching pi_M and kappa of
# LEYE_recruitment_two_era.rds reproduces Beta(1.8410, 32.5058) and
# Beta(7.6338, 8.6569) exactly -- so every setting below is built the same way
# the main fit was, not by a parallel formula that could drift from it.
beta_mom <- function(x) {
  mu <- mean(x)
  k <- mu * (1 - mu) / var(x) - 1
  stopifnot(is.finite(k), k > 0)
  c(a = mu * k, b = (1 - mu) * k)
}

# Both recruitment parameters re-derived from the fit with 1997 removed.
drop97_prior <- function() {
  f <- "LEYE_recruitment_two_era_drop1997.rds"
  if (!file.exists(f)) {
    stop("Missing ", f, ". Build it first with:\n",
         "  Rscript scripts/LEYE_recruitment_drop_cohort.R 1997")
  }
  te <- do.call(rbind, readRDS(f))
  p <- beta_mom(te[, "pi_M"])
  k <- beta_mom(te[, "kappa"])
  list(a_pi = p[["a"]], b_pi = p[["b"]], a_kap = k[["a"]], b_kap = k[["b"]])
}

# Beta fit to the cohort-averaged modern pi.
cohortmean_prior <- function() {
  te <- do.call(rbind, readRDS("LEYE_recruitment_two_era.rds"))
  pic <- grep("^pi_c\\[", colnames(te), value = TRUE)
  N_COH_90 <- 5 # 1995-1999; modern cohorts follow (see `modern` in that module)
  modern <- pic[(N_COH_90 + 1):length(pic)]
  stopifnot(length(modern) > 0)
  p <- beta_mom(rowMeans(te[, modern, drop = FALSE]))
  list(a_pi = p[["a"]], b_pi = p[["b"]])
}

SETTINGS <- list(
  cohortmean = list(
    out = "LEYE_IPM_samples_cohortmean.rds",
    label = "cohort-averaged pi",
    # kappa unchanged: it is a scalar with no cohort random effect, so it has no
    # spread across cohorts and no typical-vs-average gap to correct.
    make = cohortmean_prior
  ),
  drop97 = list(
    out = "LEYE_IPM_samples_drop1997.rds",
    label = "1997 cohort excluded, pi and kappa",
    make = drop97_prior
  ),
  vague = list(
    out = "LEYE_IPM_samples_vaguepriors.rds",
    label = "vague Beta(1,1) on both",
    make = function() list(a_pi = 1, b_pi = 1, a_kap = 1, b_kap = 1)
  )
)

## ---------------------------------------------------------------------------
## Re-stamp without refitting
## ---------------------------------------------------------------------------
# stamp_ipm_cache() used to default to the global `constants`, so companion fits
# were saved carrying the DEFAULT recruitment priors in their provenance -- the
# samples were right, the recorded priors were not, and the report printed them.
# This rewrites the stamp from each setting's real constants, keeping the
# original fit timestamp. Samples are untouched.
restamp_all <- function() {
  for (key in names(SETTINGS)) {
    s <- SETTINGS[[key]]
    if (!file.exists(s$out)) {
      cat(sprintf("  %-12s %s not found, skipped\n", key, s$out))
      next
    }
    sv <- readRDS(s$out)
    old <- attr(sv, "ipm_provenance")
    cs <- constants
    pr <- s$make()
    for (nm in names(pr)) cs[[nm]] <- pr[[nm]]
    new <- ipm_provenance(consts = cs)
    if (!is.null(old$fitted)) new$fitted <- old$fitted
    attr(sv, "ipm_provenance") <- new
    saveRDS(sv, s$out)
    changed <- is.null(old$priors) ||
      !isTRUE(all.equal(unname(old$priors), unname(new$priors)))
    cat(sprintf("  %-12s pi ~ Beta(%.4f, %.4f)  kappa ~ Beta(%.4f, %.4f)%s\n",
                key, cs$a_pi, cs$b_pi, cs$a_kap, cs$b_kap,
                if (changed) "   <- corrected" else "   (already correct)"))
  }
}

args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "restamp")) {
  cat("\n=== re-stamping companion provenance (no refit) ===\n")
  restamp_all()
  quit(save = "no")
}
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
# Overridable from the environment so a fit that needs a longer chain does not
# need a forked script. `logB0` is weakly identified against the population
# level in EVERY fit -- cor(logB0, mean log N_ad) = -0.85, ESS ~90-120 out of
# 21000 draws -- so its Rhat is noisy at the default length and can land
# anywhere from 1.01 to 1.19 on luck alone. The posterior is unchanged by a
# longer run, so a companion fitted at a different chain length stays
# comparable to the main fit; only the Monte Carlo error differs.
#   LEYE_NITER=300000 LEYE_NBURNIN=60000 LEYE_THIN=25 Rscript ... vague
niter <- as.integer(Sys.getenv("LEYE_NITER", "100000"))
nburnin <- as.integer(Sys.getenv("LEYE_NBURNIN", "30000"))
thin <- as.integer(Sys.getenv("LEYE_THIN", "10"))
stopifnot(niter > nburnin, thin >= 1, (niter - nburnin) / thin >= 1000)

run_one <- function(key) {
  s <- SETTINGS[[key]]
  cs <- constants
  pr <- s$make() # once: each call re-reads and rbinds a whole posterior
  for (nm in names(pr)) cs[[nm]] <- pr[[nm]]

  cat(sprintf(
    "\n=== %s (%s) ===\n  pi_rec ~ Beta(%.4f, %.4f)   kappa ~ Beta(%.4f, %.4f)\n",
    key, s$label, cs$a_pi, cs$b_pi, cs$a_kap, cs$b_kap
  ))
  cat(sprintf("  main fit uses Beta(%.4f, %.4f) and Beta(%.4f, %.4f)\n",
              constants$a_pi, constants$b_pi, constants$a_kap, constants$b_kap))
  cat(sprintf("  %d chains x %d iter, %d burnin, thin %d -> %d draws/chain\n",
              nchains, niter, nburnin, thin, (niter - nburnin) %/% thin))

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
  saveRDS(stamp_ipm_cache(samples, cs), s$out) # cs, not the global constants
  cat(sprintf(
    "  done in %.1f min -> %s\n",
    as.numeric(difftime(Sys.time(), t0, units = "mins")), s$out
  ))
}

for (k in which_run) run_one(k)
cat("\nAll requested sensitivity fits complete.\n")
