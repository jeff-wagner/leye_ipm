# =============================================================================
# LEYE local recruitment: two-era model
#
# Estimates the probability that a fledged chick survives AND returns to breed
# locally, separately for the 1990s and the modern era, and tests whether it
# changed. Supersedes scripts/LEYE_juv90_module.R, which fitted the 1990s alone.
#
# DATA
#   1990s chicks   139 banded 1995-1998, 6 occasions (1995-2000)
#   1990s adults   178 contemporaneous adults from the main encounter history
#   modern chicks  133 banded 2018-2025, 9 occasions (2018-2026)
#
# WHY THIS WORKS AT ALL
# Recruitment and detection are confounded in chick data alone: only their
# product is identified, so 5 modern resights out of 133 is consistent with
# anything from ~5% to ~18% recruitment depending on detection. Birds in BOTH
# eras carried colour marks read by scope, so detection is modelled as shared,
# and the 178 contemporaneous 1990s adults pin it. That is what turns the
# modern chick data from a bound into an estimate.
#
#   logit(p[t]) = mu_p + eps_p[t]        one detection process, all years
#
# Detection is partially pooled rather than forced equal: era or year effects
# can still show up, but modern years -- which have almost no information of
# their own -- are shrunk toward the level the 1990s adults establish.
#
# RECRUITMENT
#   logit(pi[c]) = mu_pi + b_era * modern[c] + eps[c]
#     b_era  the quantity of interest: change in recruitment between eras
#     eps[c] cohort random effect, absorbing between-year variation and any
#            residual error in cohort totals
#   psi1 = pi * kappa, psi2 = pi * (1 - kappa)   split between ages 1 and 2
#
# kappa is shared across eras: the age-at-first-return patterns are nearly
# identical (1990s 6/6/1 at ages 1/2/3, modern 2/2/1), and neither era alone
# would support a separate estimate.
#
# ADULT SURVIVAL
# Estimated for the 1990s from the contemporaneous adults. There are no modern
# adults in this module -- they are flagged rather than colour-marked, so they
# do not share the detection process -- and modern adult survival therefore
# carries an informative prior from the IPM posterior (phi_ad ~ 0.63-0.70).
# =============================================================================

library(nimble)
library(coda)
set.seed(20260807)

if (!exists("SKIP_RECRUIT_RUN")) SKIP_RECRUIT_RUN <- FALSE

## ===========================================================================
## SECTION 1 -- DATA
## ===========================================================================
Y90 <- 1995:2000
YM <- 2018:2026
N90 <- length(Y90)
NM <- length(YM)

# ---- 1990s chicks ----
j90 <- read.csv("data/LEYE_juv_EH_1995_2000.csv", check.names = FALSE)
C90 <- as.matrix(j90[, as.character(Y90)])
b90 <- match(j90$`Year banded`, Y90)
last90 <- apply(C90, 1, function(r) max(which(!is.na(r)))) # NA = not monitored
C90[is.na(C90)] <- 0L
storage.mode(C90) <- "integer"
keep <- b90 < last90
C90 <- C90[keep, , drop = FALSE]
coh90 <- as.integer(factor(j90$`Year banded`[keep]))
b90 <- b90[keep]
last90 <- last90[keep]

# ---- 1990s adults ----
eh <- read.csv("data/LEYE_95_26_EH.csv", check.names = FALSE)
A90 <- as.matrix(eh[, c("yr_95", "yr_96", "yr_97", "yr_98", "yr_99", "yr_00")])
A90[is.na(A90)] <- 0L
storage.mode(A90) <- "integer"
A90 <- A90[rowSums(A90 == 1) >= 1, , drop = FALSE]
fA90 <- max.col(A90 == 1, ties.method = "first")
A90 <- A90[fA90 < N90, , drop = FALSE]
fA90 <- fA90[fA90 < N90]

# ---- modern chicks ----
cm <- read.csv("data/LEYE_chick_EH_2018_2026.csv",
  check.names = FALSE, colClasses = c(band_number = "character")
)
CM <- as.matrix(cm[, as.character(YM)])
bM <- match(cm$year_banded, YM)
CM[is.na(CM)] <- 0L
storage.mode(CM) <- "integer"
keep <- bM < NM
CM <- CM[keep, , drop = FALSE]
cohM <- as.integer(factor(cm$year_banded[keep]))
bM <- bM[keep]

nC90 <- nrow(C90)
nA90 <- nrow(A90)
nCM <- nrow(CM)
nCoh90 <- max(coh90)
nCohM <- max(cohM)
# cohort index runs across both eras; `modern` flags which era each belongs to
cohM_off <- cohM + nCoh90
nCohAll <- nCoh90 + nCohM
modern <- c(rep(0, nCoh90), rep(1, nCohM))

cat(sprintf(
  "1990s: %d chicks (%d cohorts) + %d adults | modern: %d chicks (%d cohorts)\n",
  nC90, nCoh90, nA90, nCM, nCohM
))

## ===========================================================================
## SECTION 2 -- LIKELIHOOD (validated in scripts/LEYE_juv90_module.R)
## ===========================================================================
SKIP_RUN_JUV <- TRUE
suppressMessages(source("scripts/LEYE_juv90_module.R"))
# provides: chiFrom, cjsCore, dChick/rChick, dAdult/rAdult, all registered and
# checked against brute-force enumeration AND for summing to 1 over all
# possible encounter histories.

## ===========================================================================
## SECTION 3 -- MODEL
## ===========================================================================
recruitCode <- nimbleCode({
  # ---- detection: one process, partially pooled across all years -----------
  mu_p ~ dnorm(0, sd = 1.5)
  sd_p ~ dexp(2)
  for (t in 1:nYrAll) {
    eps_p[t] ~ dnorm(0, sd = sd_p)
    logit(p_all[t]) <- mu_p + eps_p[t]
  }
  for (t in 1:n90) {
    p90[t] <- p_all[t]
  }
  for (t in 1:nM) {
    pM[t] <- p_all[n90 + t]
  }

  # ---- adult survival ------------------------------------------------------
  mu_phi90 ~ dnorm(0, sd = 1.5)
  # no modern adults share this detection process, so their survival comes
  # from the IPM posterior rather than from data in this module
  mu_phiM ~ dnorm(phiM_mu, sd = phiM_sd) # default logit(0.68) ~ 0.75
  for (t in 1:(n90 - 1)) {
    phi90[t] <- ilogit(mu_phi90)
  }
  for (t in 1:(nM - 1)) {
    phiM[t] <- ilogit(mu_phiM)
  }
  phi_ad_90 <- ilogit(mu_phi90)
  phi_ad_M <- ilogit(mu_phiM)

  # ---- recruitment ---------------------------------------------------------
  mu_pi ~ dnorm(0, sd = 1.5)
  b_era ~ dnorm(0, sd = 1.5) # modern vs 1990s, logit scale
  sd_pi ~ dexp(2)
  kappa ~ dbeta(1, 1)
  for (c in 1:nCohAll) {
    eps_pi[c] ~ dnorm(0, sd = sd_pi)
    logit(pi_c[c]) <- mu_pi + b_era * modern[c] + eps_pi[c]
    psi1_c[c] <- pi_c[c] * kappa
    psi2_c[c] <- pi_c[c] * (1 - kappa)
  }
  pi_90 <- ilogit(mu_pi)
  pi_M <- ilogit(mu_pi + b_era)
  era_ratio <- pi_M / pi_90

  # ---- likelihoods ---------------------------------------------------------
  for (i in 1:nC90) {
    y90[i, 1:n90] ~ dChick(
      psi1_c[coh90[i]], psi2_c[coh90[i]],
      phi90[1:(n90 - 1)], p90[1:n90], band90[i], last90[i]
    )
  }
  for (j in 1:nA90) {
    a90[j, 1:n90] ~ dAdult(phi90[1:(n90 - 1)], p90[1:n90], fa90[j], n90)
  }
  for (i in 1:nCM) {
    yM[i, 1:nM] ~ dChick(
      psi1_c[cohMoff[i]], psi2_c[cohMoff[i]],
      phiM[1:(nM - 1)], pM[1:nM], bandM[i], nM
    )
  }
})

constants <- list(
  n90 = N90, nM = NM, nYrAll = N90 + NM,
  phiM_mu = 0.75, phiM_sd = 0.4,
  nC90 = nC90, nA90 = nA90, nCM = nCM, nCohAll = nCohAll,
  coh90 = coh90, cohMoff = cohM_off, modern = modern,
  band90 = b90, last90 = last90, fa90 = fA90, bandM = bM
)
data <- list(y90 = C90, a90 = A90, yM = CM)
inits <- function() {
  list(
    mu_p = qlogis(0.7), sd_p = 0.4, eps_p = rnorm(N90 + NM, 0, 0.2),
    mu_phi90 = qlogis(0.65), mu_phiM = qlogis(0.68),
    mu_pi = qlogis(0.12), b_era = 0, sd_pi = 0.5,
    eps_pi = rnorm(nCohAll, 0, 0.3), kappa = 0.5
  )
}
monitors <- c(
  "pi_90", "pi_M", "b_era", "era_ratio", "kappa",
  "phi_ad_90", "phi_ad_M", "mu_p", "sd_p", "sd_pi", "p_all", "pi_c"
)

## ===========================================================================
## SECTION 4 -- RUN
## ===========================================================================
if (!SKIP_RECRUIT_RUN) {
  Rm <- nimbleModel(recruitCode,
    constants = constants, data = data,
    inits = inits(), calculate = FALSE
  )
  cat("logProb at inits:", Rm$calculate(), "\n")
  conf <- configureMCMC(Rm, monitors = monitors, print = FALSE)
  mc <- buildMCMC(conf)
  Cm <- compileNimble(Rm)
  Cmc <- compileNimble(mc, project = Rm)
  samples <- runMCMC(Cmc,
    niter = 80000, nburnin = 30000, thin = 10, nchains = 3,
    inits = inits, samplesAsCodaMCMC = TRUE, setSeed = 1:3
  )

  m <- as.matrix(samples)
  q <- function(x) {
    sprintf("%.3f (%.3f, %.3f)", mean(x), quantile(x, .025), quantile(x, .975))
  }
  cat("\n=====================================================================\n")
  cat("LOCAL RECRUITMENT, TWO ERAS\n")
  cat("=====================================================================\n")
  kp <- c("pi_90", "pi_M", "b_era", "kappa", "phi_ad_90", "mu_p", "sd_p", "sd_pi")
  cat("max Rhat:", round(max(gelman.diag(samples[, kp], multivariate = FALSE)$psrf[, 1]), 4), "\n\n")
  cat("recruitment, 1990s      :", q(m[, "pi_90"]), "\n")
  cat("recruitment, modern     :", q(m[, "pi_M"]), "\n")
  cat("era effect (logit)      :", q(m[, "b_era"]), "\n")
  cat("ratio modern / 1990s    :", q(m[, "era_ratio"]), "\n")
  cat(sprintf("P(recruitment declined) : %.3f\n", mean(m[, "b_era"] < 0)))
  cat("\nfraction recruiting at age 1 (kappa):", q(m[, "kappa"]), "\n")
  cat("1990s adult survival    :", q(m[, "phi_ad_90"]), "\n")
  cat("detection, mean over yrs:", q(rowMeans(m[, paste0("p_all[", 1:(N90 + NM), "]")])), "\n")
  cat("  1990s years           :", q(rowMeans(m[, paste0("p_all[", 1:N90, "]")])), "\n")
  cat("  modern years          :", q(rowMeans(m[, paste0("p_all[", (N90 + 1):(N90 + NM), "]")])), "\n")

  saveRDS(samples, "LEYE_recruitment_two_era.rds")
  cat("\nWrote: LEYE_recruitment_two_era.rds\n")
}
