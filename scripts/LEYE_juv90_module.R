# =============================================================================
# LEYE 1990s JUVENILE RECRUITMENT MODULE -- NIMBLE
#
# Purpose: estimate local recruitment probabilities for chicks banded 1995-1998,
# to replace the weakly-identified `rho` (juvenile:adult survival ratio) in the
# main IPM (scripts/LEYE_IPM.R).
#
# Two datasets, sharing adult survival and year-specific detection:
#   (1) CHICKS  : 139 birds banded as chicks 1995-1998 (data/LEYE_juv_EH_1995_2000.csv)
#   (2) ADULTS  : 178 contemporaneous adults from the main EH, occasions 1995-2000
#
# Why a custom likelihood: a chick first SEEN at age 2 may have recruited at age
# 1 and gone undetected. Raw return rates confound recruitment timing with
# detection. Here recruitment age is a latent mixture, marginalised analytically:
#
#   A_i = 1    w.p. psi1[cohort]      (recruits locally at age 1)
#   A_i = 2    w.p. psi2[cohort]      (recruits locally at age 2)
#   A_i = never w.p. 1 - psi1 - psi2
#
# psi1/psi2 are COMPOSITE probabilities (survive AND return to the local
# breeding population), i.e. exactly the quantity the IPM recruitment equation
# needs. They are not annual survival rates -- do not give them survival priors.
#
# Parameterisation:
#   logit(pi[c]) = mu_pi + eps[c]     total recruitment prob, cohort random effect
#   psi1[c] = pi[c] * kappa           kappa = fraction recruiting at age 1
#   psi2[c] = pi[c] * (1 - kappa)
#
# The cohort random effect also absorbs uncertainty in the approximate cohort
# totals (37/34/20/48), which enter the likelihood the same way as genuine
# cohort heterogeneity.
#
# NA in the chick histories = not monitored (right-censoring), handled via a
# per-individual `last monitored occasion`.
# =============================================================================

library(nimble)
library(coda)
set.seed(20260804)

## ===========================================================================
## SECTION 1 -- DATA PREP
## ===========================================================================
NOCC <- 6L # occasions 1995..2000
occ_years <- 1995:2000

# ---- 1a. Chicks ------------------------------------------------------------
jv <- read.csv("data/LEYE_juv_EH_1995_2000.csv", check.names = FALSE)
ycols <- as.character(occ_years)
Ych <- as.matrix(jv[, ycols]) # may contain NA = not monitored

band_occ <- match(jv$`Year banded`, occ_years)
stopifnot(!anyNA(band_occ))

# last monitored occasion = last non-NA entry (NAs run to end of study)
last_occ <- apply(Ych, 1, function(r) max(which(!is.na(r))))
stopifnot(all(last_occ >= band_occ))

# after censoring bookkeeping, NAs never get read; zero them so the matrix is numeric
Ych[is.na(Ych)] <- 0
storage.mode(Ych) <- "integer"

cohort <- as.integer(factor(jv$`Year banded`))
cohort_years <- sort(unique(jv$`Year banded`))
nCoh <- length(cohort_years)
nChick <- nrow(Ych)

# chicks whose banding occasion is the final occasion contribute nothing
keep_ch <- band_occ < last_occ
Ych <- Ych[keep_ch, , drop = FALSE]
band_occ <- band_occ[keep_ch]
last_occ <- last_occ[keep_ch]
cohort <- cohort[keep_ch]
nChick <- nrow(Ych)

# ---- 1b. Contemporaneous adults from the main EH ---------------------------
eh <- read.csv("data/LEYE_95_26_EH.csv", check.names = FALSE)
cv <- read.csv("data/LEYE_95_26_covs.csv", check.names = FALSE)
stopifnot(identical(eh$Band_Number, cv$Band_Number)) # lockstep required

ad_cols <- c("yr_95", "yr_96", "yr_97", "yr_98", "yr_99", "yr_00")
Yad <- as.matrix(eh[, ad_cols])
Yad[is.na(Yad)] <- 0
storage.mode(Yad) <- "integer"
Yad <- Yad[rowSums(Yad == 1) >= 1, , drop = FALSE]
fad <- max.col(Yad == 1, ties.method = "first")
Yad <- Yad[fad < NOCC, , drop = FALSE] # need >=1 interval at risk
fad <- fad[fad < NOCC]
nAd <- nrow(Yad)

cat(sprintf(
  "Chicks: %d (cohorts %s)   Adults: %d   Occasions: %d\n",
  nChick,
  paste(cohort_years, collapse = "/"),
  nAd,
  NOCC
))

## ===========================================================================
## SECTION 2 -- LIKELIHOOD BUILDING BLOCKS
## ===========================================================================

# chi[t0] = P(no detection in (t0, last] | alive at t0)
chiFrom <- nimbleFunction(
  run = function(
    t0 = double(0),
    last = double(0),
    phi = double(1),
    p = double(1)
  ) {
    returnType(double(0))
    ch <- 1
    if (last > t0) {
      for (k in 1:(last - t0)) {
        t <- last - k
        ch <- (1 - phi[t]) + phi[t] * (1 - p[t + 1]) * ch
      }
    }
    return(ch)
  }
)

# CJS likelihood over occasions [e, last], given alive & available at e.
# entryDet = 1 -> the detection outcome AT e is informative (chicks: availability
#                 is governed by recruitment, not by having been captured)
# entryDet = 0 -> condition on release at e, no detection term at e (adults)
cjsCore <- nimbleFunction(
  run = function(
    x = double(1),
    e = double(0),
    last = double(0),
    phi = double(1),
    p = double(1),
    entryDet = double(0)
  ) {
    returnType(double(0))
    L <- 0
    for (t in e:last) {
      if (x[t] == 1) L <- t
    }
    if (L == 0) {
      # never detected on [e, last]
      base <- chiFrom(e, last, phi, p)
      if (entryDet == 1) {
        return((1 - p[e]) * base)
      }
      return(base)
    }
    val <- 1
    if (entryDet == 1) {
      if (x[e] == 1) {
        val <- val * p[e]
      } else {
        val <- val * (1 - p[e])
      }
    }
    if (L > e) {
      for (t in (e + 1):L) {
        if (x[t] == 1) {
          val <- val * p[t]
        } else {
          val <- val * (1 - p[t])
        }
      }
      for (t in e:(L - 1)) {
        val <- val * phi[t]
      }
    }
    return(val * chiFrom(L, last, phi, p))
  }
)

# --- chick distribution: mixture over recruitment age -----------------------
dChick <- nimbleFunction(
  run = function(
    x = double(1),
    psi1 = double(0),
    psi2 = double(0),
    phi = double(1),
    p = double(1),
    band = double(0),
    last = double(0),
    log = integer(0, default = 0)
  ) {
    returnType(double(0))
    lik <- 0
    # component: never recruits locally -> cannot be detected after banding
    anydet <- 0
    for (t in (band + 1):last) {
      if (x[t] == 1) anydet <- 1
    }
    if (anydet == 0) lik <- lik + (1 - psi1 - psi2)
    # components: recruits at age 1 or 2
    for (a in 1:2) {
      e <- band + a
      psia <- psi1
      if (a == 2) psia <- psi2
      if (e > last) {
        # recruits after the observation window closes -> unobservable here
        if (anydet == 0) lik <- lik + psia
      } else {
        # a bird recruiting at age `a` is NOT in the population before then,
        # so it cannot have been detected on (band, e). Without this guard the
        # age-2 component ignores occasion band+1 and double-counts it.
        predet <- 0
        if (e > band + 1) {
          for (t in (band + 1):(e - 1)) {
            if (x[t] == 1) predet <- 1
          }
        }
        if (predet == 0) lik <- lik + psia * cjsCore(x, e, last, phi, p, 1)
      }
    }
    if (log) return(log(lik))
    return(lik)
  }
)
rChick <- nimbleFunction(
  run = function(
    n = integer(0),
    psi1 = double(0),
    psi2 = double(0),
    phi = double(1),
    p = double(1),
    band = double(0),
    last = double(0)
  ) {
    returnType(double(1))
    return(numeric(length = 6)) # not used; required for registration
  }
)

# --- adult distribution: standard CJS, conditioned on release ---------------
dAdult <- nimbleFunction(
  run = function(
    x = double(1),
    phi = double(1),
    p = double(1),
    f = double(0),
    last = double(0),
    log = integer(0, default = 0)
  ) {
    returnType(double(0))
    lik <- cjsCore(x, f, last, phi, p, 0)
    if (log) return(log(lik))
    return(lik)
  }
)
rAdult <- nimbleFunction(
  run = function(
    n = integer(0),
    phi = double(1),
    p = double(1),
    f = double(0),
    last = double(0)
  ) {
    returnType(double(1))
    return(numeric(length = 6))
  }
)

registerDistributions(list(
  dChick = list(
    BUGSdist = "dChick(psi1, psi2, phi, p, band, last)",
    types = c(
      "value = double(1)",
      "phi = double(1)",
      "p = double(1)"
    ),
    discrete = TRUE
  ),
  dAdult = list(
    BUGSdist = "dAdult(phi, p, f, last)",
    types = c(
      "value = double(1)",
      "phi = double(1)",
      "p = double(1)"
    ),
    discrete = TRUE
  )
))

## ===========================================================================
## SECTION 3 -- MODEL
## ===========================================================================
juvCode <- nimbleCode({
  # --- 1990s adult survival (constant across the 5 intervals) ---
  mu_phi_a ~ dnorm(0, sd = 1.5)
  for (t in 1:(nocc - 1)) {
    phi_a[t] <- ilogit(mu_phi_a)
  }
  phi_ad90 <- ilogit(mu_phi_a)

  # --- year-specific detection (occasion 1 never used) ---
  p[1] <- 0
  for (t in 2:nocc) {
    p[t] ~ dbeta(1, 1)
  }

  # --- local recruitment: total probability with cohort random effect ---
  mu_pi ~ dnorm(0, sd = 1.5)
  sd_pi ~ dexp(2)
  for (c in 1:nCoh) {
    eps_pi[c] ~ dnorm(0, sd = sd_pi)
    logit(pi_c[c]) <- mu_pi + eps_pi[c]
    psi1_c[c] <- pi_c[c] * kappa
    psi2_c[c] <- pi_c[c] * (1 - kappa)
  }
  kappa ~ dbeta(1, 1) # fraction of recruits entering at age 1

  # population-level (cohort-average) quantities for the IPM
  pi_bar <- ilogit(mu_pi)
  psi1_bar <- pi_bar * kappa
  psi2_bar <- pi_bar * (1 - kappa)
  # logit contrast vs 1990s adult survival (comparable to literature)
  delta <- logit(pi_bar) - logit(phi_ad90)

  # --- likelihoods ---
  for (i in 1:nChick) {
    ych[i, 1:nocc] ~ dChick(
      psi1_c[coh[i]],
      psi2_c[coh[i]],
      phi_a[1:(nocc - 1)],
      p[1:nocc],
      band[i],
      lastc[i]
    )
  }
  for (j in 1:nAd) {
    yad[j, 1:nocc] ~ dAdult(phi_a[1:(nocc - 1)], p[1:nocc], fad[j], nocc)
  }
})

constants <- list(
  nocc = NOCC,
  nChick = nChick,
  nAd = nAd,
  nCoh = nCoh,
  coh = cohort,
  band = band_occ,
  lastc = last_occ,
  fad = fad
)
data <- list(ych = Ych, yad = Yad)
inits <- function() {
  list(
    mu_phi_a = qlogis(0.65),
    p = c(NA, runif(NOCC - 1, 0.3, 0.7)),
    mu_pi = qlogis(0.15),
    sd_pi = 0.5,
    eps_pi = rnorm(nCoh, 0, 0.3),
    kappa = 0.5
  )
}
monitors <- c(
  "phi_ad90",
  "p",
  "pi_bar",
  "psi1_bar",
  "psi2_bar",
  "kappa",
  "delta",
  "mu_pi",
  "sd_pi",
  "pi_c"
)

## ===========================================================================
## SECTION 4 -- RUN
## ===========================================================================
if (!exists("SKIP_RUN_JUV")) SKIP_RUN_JUV <- FALSE
if (!SKIP_RUN_JUV) {
  Rmodel <- nimbleModel(
    juvCode,
    constants = constants,
    data = data,
    inits = inits(),
    calculate = FALSE
  )
  cat("logProb at inits:", Rmodel$calculate(), "\n")
  conf <- configureMCMC(Rmodel, monitors = monitors)
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

  samples <- runMCMC(
    Cmcmc,
    niter = 60000,
    nburnin = 20000,
    thin = 10,
    nchains = 3,
    inits = inits,
    samplesAsCodaMCMC = TRUE,
    setSeed = 1:3
  )

  library(MCMCvis)
  print(MCMCsummary(samples, round = 3))

  mat <- as.matrix(samples)
  q <- function(x) {
    c(
      mean = mean(x),
      lwr = quantile(x, .025, names = FALSE),
      upr = quantile(x, .975, names = FALSE)
    )
  }
  cat("\n--- Quantities for the IPM ---\n")
  out <- rbind(
    psi1 = q(mat[, "psi1_bar"]),
    psi2 = q(mat[, "psi2_bar"]),
    pi_total = q(mat[, "pi_bar"]),
    kappa = q(mat[, "kappa"]),
    phi_ad_1990s = q(mat[, "phi_ad90"]),
    delta_logit = q(mat[, "delta"])
  )
  print(round(out, 3))
  saveRDS(samples, "LEYE_juv90_samples.rds")
}
