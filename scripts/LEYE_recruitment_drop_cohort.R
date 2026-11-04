# =============================================================================
# LEYE local recruitment: leave-one-cohort-out sensitivity
#
# Refits scripts/LEYE_recruitment_two_era.R with one or more natal-year cohorts
# deleted, to measure how much of the era contrast rests on a single cohort.
#
# WHY
# The 1997 cohort returns pi_c = 0.584 (0.252-0.884) against <=0.12 for every
# other cohort in either era, on 11 chicks banded and 7 resighted. It is also
# the cohort that sets sd_pi = 1.32 on the logit scale, and sd_pi is shared
# across eras -- so one cohort inflates the uncertainty on all thirteen.
#
# This is deliberately NOT a brood random effect. The seven 1997 returns carry
# seven distinct brood/site prefixes (Connor's, Turnagain, Nordic, Crane,
# Detox, Fancy, Horsefly), so none of them are siblings and brood pooling
# would not touch the anomaly. Brood is also unassignable for most of the
# 1990s: 106 of the 122 chicks are reconstructed never-resighted rows carrying
# nothing but natal year, and recorded brood IDs survive for 1996 only.
# Dropping the cohort tests its influence directly instead.
#
# WHAT CHANGES
# Only the chick rows in the dropped cohorts, and the cohort index they share.
# Priors, detection, adult survival, kappa and both likelihoods are inherited
# from the parent script unchanged -- the 1990s adults stay in, so detection is
# still pinned by them.
#
# USAGE  (from the PROJECT ROOT)
#   Rscript scripts/LEYE_recruitment_drop_cohort.R          # drops 1997
#   Rscript scripts/LEYE_recruitment_drop_cohort.R 1997 1996
#
# WRITES  LEYE_recruitment_two_era_drop<years>.rds
# =============================================================================

DROP <- as.integer(commandArgs(trailingOnly = TRUE))
if (!length(DROP) || anyNA(DROP)) DROP <- 1997L

# Inherit data, model code and priors from the parent with its own run skipped.
SKIP_RECRUIT_RUN <- TRUE
source("scripts/LEYE_recruitment_two_era.R")

suppressMessages(library(coda))

stopifnot(all(DROP %in% coh_years_all))
tag <- paste(DROP, collapse = "_")
cat(sprintf("\n==== dropping cohort(s): %s ====\n", paste(DROP, collapse = ", ")))

## ---- rebuild the cohort index without the dropped years ---------------------
keep_coh <- !(coh_years_all %in% DROP)     # over 1:nCohAll
new_idx <- cumsum(keep_coh)                 # old cohort index -> new, where kept

d90 <- which(coh90_years %in% DROP)         # 1990s cohort positions dropped
dM <- which(cohM_years %in% DROP) + nCoh90  # modern, on the combined scale

keep90 <- !(coh90 %in% d90)
keepM <- !(cohM_off %in% dM)

C90d <- C90[keep90, , drop = FALSE]
CMd <- CM[keepM, , drop = FALSE]
coh90d <- new_idx[coh90[keep90]]
cohMd <- new_idx[cohM_off[keepM]]
modernd <- modern[keep_coh]
nCohAlld <- sum(keep_coh)

cat(sprintf(
  "chicks: 1990s %d -> %d | modern %d -> %d | cohorts %d -> %d\n",
  nC90, nrow(C90d), nCM, nrow(CMd), nCohAll, nCohAlld
))
stopifnot(
  max(coh90d, cohMd) == nCohAlld,
  length(unique(c(coh90d, cohMd))) == nCohAlld, # no cohort left empty
  sum(modernd) == nCohM - length(dM), nA90 == nrow(A90)
)

constants_d <- constants
constants_d[c("nC90", "nCM", "nCohAll", "coh90", "cohMoff", "modern",
  "band90", "last90", "bandM")] <- list(
  nrow(C90d), nrow(CMd), nCohAlld, coh90d, cohMd, modernd,
  b90[keep90], last90[keep90], bM[keepM]
)
data_d <- list(y90 = C90d, a90 = A90, yM = CMd)
inits_d <- function() {
  x <- inits()
  x$eps_pi <- rnorm(nCohAlld, 0, 0.3)
  x
}

## ---- fit --------------------------------------------------------------------
Rm <- nimbleModel(recruitCode,
  constants = constants_d, data = data_d,
  inits = inits_d(), calculate = FALSE
)
cat("logProb at inits:", Rm$calculate(), "\n")
stopifnot(is.finite(Rm$calculate()))
conf <- configureMCMC(Rm, monitors = monitors, print = FALSE)
mc <- buildMCMC(conf)
Cm <- compileNimble(Rm)
Cmc <- compileNimble(mc, project = Rm)
samples <- runMCMC(Cmc,
  niter = 80000, nburnin = 30000, thin = 10, nchains = 3,
  inits = inits_d, samplesAsCodaMCMC = TRUE, setSeed = 1:3
)

m <- as.matrix(samples)
q <- function(x) {
  sprintf("%.3f (%.3f, %.3f)", mean(x), quantile(x, .025), quantile(x, .975))
}
kp <- c("pi_90", "pi_M", "b_era", "kappa", "phi_ad_90", "mu_p", "sd_p", "sd_pi")

cat("\n=====================================================================\n")
cat(sprintf("LOCAL RECRUITMENT, TWO ERAS -- cohort(s) %s removed\n", tag))
cat("=====================================================================\n")
cat("max Rhat:", round(max(gelman.diag(samples[, kp],
  multivariate = FALSE
)$psrf[, 1]), 4), "\n\n")
cat("recruitment, 1990s      :", q(m[, "pi_90"]), "\n")
cat("recruitment, modern     :", q(m[, "pi_M"]), "\n")
cat("era effect (logit)      :", q(m[, "b_era"]), "\n")
cat("ratio modern / 1990s    :", q(m[, "era_ratio"]), "\n")
cat(sprintf("P(recruitment declined) : %.3f\n", mean(m[, "b_era"] < 0)))
cat("cohort SD (logit)       :", q(m[, "sd_pi"]), "\n")
cat("fraction at age 1 (kap) :", q(m[, "kappa"]), "\n")
cat("1990s adult survival    :", q(m[, "phi_ad_90"]), "\n")

pic <- grep("pi_c[", colnames(m), fixed = TRUE, value = TRUE)
cat("\ncohort-level pi_c:\n")
lab <- coh_years_all[keep_coh]
for (i in seq_along(pic)) cat(sprintf("  %-6s %s\n", lab[i], q(m[, pic[i]])))

## ---- side by side with the full fit ----------------------------------------
full_path <- "LEYE_recruitment_two_era.rds"
if (file.exists(full_path)) {
  fm <- as.matrix(readRDS(full_path))
  cat("\n---- full fit vs cohort(s) removed ----------------------------------\n")
  cat(sprintf("%-24s %-26s %-26s\n", "", "all cohorts", paste("drop", tag)))
  for (p in c("pi_90", "pi_M", "b_era", "era_ratio", "sd_pi", "kappa")) {
    cat(sprintf("%-24s %-26s %-26s\n", p, q(fm[, p]), q(m[, p])))
  }
  cat(sprintf(
    "%-24s %-26.3f %-26.3f\n", "P(declined)",
    mean(fm[, "b_era"] < 0), mean(m[, "b_era"] < 0)
  ))
}

out <- sprintf("LEYE_recruitment_two_era_drop%s.rds", tag)
saveRDS(samples, out)
cat("\nWrote:", out, "\n")
