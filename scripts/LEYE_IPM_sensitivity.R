# =============================================================================
# IPM sensitivity runs
#
# Refits the IPM with ONE thing changed and everything else -- structure,
# sampler settings, the other two data streams -- exactly as scripts/LEYE_IPM.R
# currently defines it, so these fits stay comparable to the main one by
# construction rather than by anyone remembering to keep them in step.
#
# Two families of setting, and they answer different questions:
#
#   RECRUITMENT PRIORS  cohortmean, drop97, vague
#     Only the Beta hyperparameters move. Data untouched.
#
#   COUNT STREAM        full, balanced9
#     Only which survey wetlands enter the count likelihood moves. Priors
#     untouched. These exist because the count stream is drawn from a wider
#     area, and a less stable panel, than the two demographic streams:
#
#       - The 32 surveyed wetlands span 67 km; five of them (Palmer Hay Flats
#         x2, Eklutna, Big Lake, Lynn Fuller's Ponds) sit 37-67 km from
#         Anchorage and contribute 22.6% of all birds counted, while every
#         marked bird and every monitored nest comes from the Anchorage bowl.
#         SINCE 2026-11-04 THE MAIN MODEL DROPS THOSE FIVE -- see the scope
#         block in LEYE_IPM.R -- so `full` is now the sensitivity in the other
#         direction: what including the wider area does.
#       - The panel GREW from 13 wetlands in 2014 to 32 by 2021, and only 9
#         were surveyed in all 13 years. The wetlands added late are the rich
#         ones: correlation between a site's mean count and its mean survey
#         year is +0.50, and the two richest sites of all (Coastal Refuge at
#         7.8 birds/survey, Oceanview Ghost Forest at 5.2) were surveyed in
#         only 7 and 8 years, centred on 2023. `balanced9` keeps the nine
#         wetlands present throughout.
#
#     The model already gives every wetland its own level offset -- `site_re`,
#     with sd_site about 1.5 on the log scale -- so it is not level that is at
#     issue. What it assumes is that the wetlands it keeps share ONE temporal
#     trajectory, `N_tot[t]`.
#
#     Compare lambda and the SHAPE of N_tot, never logB0 or absolute level:
#     logB0 is weakly identified against the population level in every fit
#     (see the chain-length note below). scripts/compare_count_subsets.R does
#     this and also checks whether site_re tracks observed richness.
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
#   full        all 32 wetlands, i.e. the count stream as it was before the
#               main model was restricted to the bowl. 935 records.
#   balanced9   only the nine wetlands surveyed in all 13 years, so panel
#               composition is constant. 333 records. NOT nested in the main
#               fit -- it keeps Eklutna and Reflections Lake, which are outside
#               the bowl but were surveyed throughout -- so the two cuts stay
#               independent: the main model's restriction isolates geography,
#               balanced9 isolates whether the changing panel manufactures
#               trend.
#
# USAGE
#   Rscript scripts/LEYE_IPM_sensitivity.R              # every setting (slow)
#   Rscript scripts/LEYE_IPM_sensitivity.R cohortmean   # just one
#   Rscript scripts/LEYE_IPM_sensitivity.R core27 balanced9   # several
#   Rscript scripts/LEYE_IPM_sensitivity.R dryrun       # validate settings, no fit
#   Rscript scripts/LEYE_IPM_sensitivity.R restamp      # fix provenance, no refit
#
# Each setting uses three cores, so on a many-core machine the count-stream
# fits are best launched as separate concurrent processes rather than as one
# invocation that works through them in series.
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

## ---------------------------------------------------------------------------
## Count-stream subsets
## ---------------------------------------------------------------------------
# A setting declares WHICH WETLANDS it keeps and nothing else; all the index
# bookkeeping happens once, here, where it can be asserted. Spelling the
# overrides out per setting would mean five chances per setting to mis-index
# `cyr` against N_tot and get a fit that runs perfectly and means nothing.
#
# Three things are deliberately NOT re-derived from the subset:
#   * `cyr` is matched against the FULL count_years, so year k still means
#     N_tot[k]. Re-factoring would silently renumber years if a subset ever
#     lost one -- asserted against below rather than trusted.
#   * `nWeek` and weekFE stay at 3; asserted.
#   * `log_et` / `log_ed` keep their ORIGINAL centring, on the full 935
#     records. Re-centring per subset would move b_et, b_ed and logB0 onto a
#     different scale and make the fits incomparable in exactly the terms we
#     want to compare them in.
# NB: everything here indexes `dc_all`, the unrestricted 935-record file, not
# `dc`. LEYE_IPM.R already restricts `dc` to the Anchorage bowl for the main
# fit, so subsetting `dc` could only ever shrink the count stream further --
# the `full` setting, which has to ADD the five distant wetlands back, would be
# impossible to express. `log_et_all` / `log_ed_all` are the matching
# all-record effort vectors, centred once there.
count_subset <- function(keep) {
  i <- which(dc_all$wetland %in% keep)
  yr <- match(dc_all$year[i], count_years)
  wk <- as.integer(factor(dc_all$week[i]))
  wet <- as.integer(factor(dc_all$wetland[i]))
  stopifnot(
    length(i) > 0,
    length(keep) == max(wet), # every requested wetland actually present
    setequal(yr, seq_len(T)), # no year lost -> N_tot indexing intact
    setequal(wk, seq_len(nWeek)) # no week lost -> weekFE indexing intact
  )
  list(
    consts = list(
      nCount = length(i), nWet = max(wet),
      cyr = yr, cwk = wk, cwet = wet
    ),
    data = list(
      count = dc_all$count[i], log_et = log_et_all[i], log_ed = log_ed_all[i]
    )
  )
}

# The nine wetlands surveyed in every year of the series, derived rather than
# listed: a hard-coded set would go stale the first time a season is added.
# Two of them (Eklutna, Reflections Lake) lie outside the Anchorage bowl, so
# this panel is NOT a subset of the main fit's -- see the setting's note.
balanced_wetlands <- function() {
  pres <- table(dc_all$wetland, dc_all$year) > 0
  as.integer(rownames(pres)[rowSums(pres) == ncol(pres)])
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
  ),
  full = list(
    out = "LEYE_IPM_samples_full.rds",
    label = "counts: all 32 wetlands, 67 km of upper Cook Inlet",
    keep = function() sort(unique(dc_all$wetland))
  ),
  balanced9 = list(
    out = "LEYE_IPM_samples_balanced9.rds",
    label = "counts: the 9 wetlands surveyed in all 13 years",
    # Not nested inside the main fit: it keeps Eklutna and Reflections Lake,
    # which are outside the bowl but were surveyed throughout. That is the
    # point -- it isolates panel COMPOSITION, where the main fit's restriction
    # isolates geography, and confounding the two would answer neither.
    keep = balanced_wetlands
  )
)
# `core27` is deliberately absent. The Anchorage bowl is now what LEYE_IPM.R
# fits, so it is the MAIN model and lives in LEYE_IPM_samples.rds. It was a
# companion until 2026-11-04; that cache was promoted rather than refitted,
# having been verified bit-identical to a run_parallel fit of the same model,
# and the duplicate was deleted so the two cannot drift.

## ---------------------------------------------------------------------------
## Re-stamp without refitting
## ---------------------------------------------------------------------------
# stamp_ipm_cache() used to default to the global `constants`, so companion fits
# were saved carrying the DEFAULT recruitment priors in their provenance -- the
# samples were right, the recorded priors were not, and the report printed them.
# This rewrites the stamp from each setting's real constants, keeping the
# original fit timestamp. Samples are untouched.
# What a setting actually fits: the global constants and data with that
# setting's overrides applied. ONE definition, used by both the runner and the
# re-stamper, so a fit and its provenance stamp cannot be built from different
# recipes -- which is the bug this section exists to clean up after.
setting_spec <- function(s) {
  cs <- constants
  dt <- data
  if (!is.null(s$make)) {
    pr <- s$make()
    for (nm in names(pr)) cs[[nm]] <- pr[[nm]]
  }
  if (!is.null(s$keep)) {
    sub <- count_subset(s$keep())
    for (nm in names(sub$consts)) cs[[nm]] <- sub$consts[[nm]]
    for (nm in names(sub$data)) dt[[nm]] <- sub$data[[nm]]
  }
  list(consts = cs, data = dt)
}

restamp_all <- function() {
  for (key in names(SETTINGS)) {
    s <- SETTINGS[[key]]
    if (!file.exists(s$out)) {
      cat(sprintf("  %-12s %s not found, skipped\n", key, s$out))
      next
    }
    sv <- readRDS(s$out)
    old <- attr(sv, "ipm_provenance")
    sp <- setting_spec(s)
    cs <- sp$consts
    new <- ipm_provenance(consts = cs, dat = sp$data)
    if (!is.null(old$fitted)) new$fitted <- old$fitted
    attr(sv, "ipm_provenance") <- new
    saveRDS(sv, s$out)
    changed <- is.null(old$priors) ||
      !isTRUE(all.equal(unname(old$priors), unname(new$priors)))
    cat(sprintf(
      "  %-12s pi ~ Beta(%.4f, %.4f)  kappa ~ Beta(%.4f, %.4f)  %d counts / %d wetlands%s\n",
      key, cs$a_pi, cs$b_pi, cs$a_kap, cs$b_kap, cs$nCount, cs$nWet,
      if (changed) "   <- corrected" else ""
    ))
  }
}

args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "restamp")) {
  cat("\n=== re-stamping companion provenance (no refit) ===\n")
  restamp_all()
  quit(save = "no")
}
# Build every setting's constants and data and print what changed, without
# fitting anything. The assertions in count_subset() are the point: an
# indexing mistake there produces a fit that runs for twenty minutes and is
# quietly meaningless, so it is worth a two-second check first.
if (identical(args, "dryrun")) {
  cat("\n=== setting dry run (no fitting) ===\n")
  cat(sprintf("  %-12s %7s %7s  %-9s %s\n",
              "setting", "counts", "wetland", "priors", "label"))
  for (key in names(SETTINGS)) {
    s <- SETTINGS[[key]]
    sp <- setting_spec(s)
    cs <- sp$consts
    same_pr <- isTRUE(all.equal(
      unname(unlist(cs[IPM_PRIOR_KEYS])), unname(unlist(constants[IPM_PRIOR_KEYS]))
    ))
    stopifnot(
      length(sp$data$count) == cs$nCount,
      length(sp$data$log_et) == cs$nCount,
      length(cs$cyr) == cs$nCount,
      setequal(cs$cyr, seq_len(T)), setequal(cs$cwk, seq_len(cs$nWeek))
    )
    cat(sprintf("  %-12s %7d %7d  %-9s %s\n", key, cs$nCount, cs$nWet,
                if (same_pr) "default" else "changed", s$label))
  }
  cat("\nAll settings build and index correctly.\n")
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
  sp <- setting_spec(s) # once: a prior setting re-reads a whole posterior
  cs <- sp$consts
  dt <- sp$data

  cat(sprintf(
    "\n=== %s (%s) ===\n  pi_rec ~ Beta(%.4f, %.4f)   kappa ~ Beta(%.4f, %.4f)\n",
    key, s$label, cs$a_pi, cs$b_pi, cs$a_kap, cs$b_kap
  ))
  cat(sprintf("  main fit uses Beta(%.4f, %.4f) and Beta(%.4f, %.4f)\n",
              constants$a_pi, constants$b_pi, constants$a_kap, constants$b_kap))
  cat(sprintf("  counts: %d records over %d wetlands (main fit %d over %d)\n",
              cs$nCount, cs$nWet, constants$nCount, constants$nWet))
  cat(sprintf("  %d chains x %d iter, %d burnin, thin %d -> %d draws/chain\n",
              nchains, niter, nburnin, thin, (niter - nburnin) %/% thin))

  # inits() closes over the GLOBAL nWet, nCount and dc$count, so a count-subset
  # fit would be handed a 32-long z_site, a 935-long od and a 935-long
  # count_rep against a model expecting fewer. Patch the three count-shaped
  # entries rather than rewriting inits() in LEYE_IPM.R, which the main fit and
  # the parallel wrapper also use.
  inits_per_chain <- lapply(seeds, function(x) {
    set.seed(x)
    i0 <- inits()
    if (!is.null(s$keep)) {
      i0$z_site <- rnorm(cs$nWet, 0, 1)
      i0$od <- rnorm(cs$nCount, 0, 0.3)
      i0$count_rep <- dt$count
    }
    i0
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
    seeds = seeds, ipmCode = ipmCode, constants = cs, data = dt,
    inits_per_chain = inits_per_chain, monitors = monitors,
    niter = niter, nburnin = nburnin, thin = thin
  )

  samples <- as.mcmc.list(chains) |>
    lapply(function(x) replace(x, is.na(x), 0)) |>
    as.mcmc.list()
  # cs and dt, never the globals: both can differ from them here.
  saveRDS(stamp_ipm_cache(samples, cs, dt), s$out)
  cat(sprintf(
    "  done in %.1f min -> %s\n",
    as.numeric(difftime(Sys.time(), t0, units = "mins")), s$out
  ))
}

for (k in which_run) run_one(k)
cat("\nAll requested sensitivity fits complete.\n")
