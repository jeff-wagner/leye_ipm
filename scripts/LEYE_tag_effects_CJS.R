# =============================================================================
# LEYE tag effects -- standalone CJS robustness analysis
#
# The IPM estimates tag effects on adult survival inside a model that also
# carries counts, productivity and a latent population process. This script
# re-estimates them from the mark-resight data ALONE, and tests two confounds
# that could manufacture an apparent tag effect. Reviewers should not have to
# accept the whole IPM to accept the tag result.
#
# Three models, all marginalising over latent alive-states (fast, no z sampling):
#
#   M0  sex + tag                      -- baseline, should reproduce the IPM
#   M1  + transience                   -- separate survival for the FIRST interval
#                                         after banding (capture effect / transients)
#   M2  + tag x time-since-banding      -- does the tag effect grow with duration?
#
# WHY THESE TWO CONFOUNDS
#
# Transience: untagged birds carry a much larger share of first-intervals
# (25.2%) than geolocator (9.0%) or GPS (14.0%) birds, because tagged cohorts
# were banded in narrow windows and followed for longer. If survival dips in
# the first interval after capture, that biases the UNTAGGED class downward and
# makes the tag effect look SMALLER than it is. M1 tests it.
#
# Duration: Pakanen et al. (2020, J Avian Biol) found shorebird survival
# declines with time carrying a device. The IPM assumes a constant tag effect.
# M2 tests a linear trend in years since banding, for tagged birds only.
#
# LIMITATION on M2: geolocators must be RECOVERED to yield data, so some birds
# stopped carrying tags before the study ended. Without per-bird recovery dates
# this is time-since-BANDING, not time-carrying. It is an attenuated proxy --
# a real duration effect would be underestimated, not manufactured.
# =============================================================================

library(nimble)
library(coda)
set.seed(20260806)

## ===========================================================================
## SECTION 1 -- DATA
## ===========================================================================
eh <- read.csv("data/LEYE_95_26_EH.csv", check.names = FALSE)
cov <- read.csv("data/LEYE_95_26_covs.csv", check.names = FALSE)
stopifnot(identical(eh$Band_Number, cov$Band_Number))

yr_of <- function(c) {
  s <- sub("yr_", "", c)
  n <- as.integer(s)
  if (n >= 93 && n <= 99) 1900 + n else if (nchar(s) == 4) n else 2000 + n
}
ycols <- grep("^yr_", names(eh), value = TRUE)
ymap <- setNames(ycols, vapply(ycols, yr_of, numeric(1)))
YEARS <- 2010:2026
M <- as.matrix(eh[, ymap[as.character(YEARS)]])
M[is.na(M)] <- 0
storage.mode(M) <- "integer"

has <- rowSums(M == 1) >= 1
M <- M[has, , drop = FALSE]
cov <- cov[has, , drop = FALSE]
f <- max.col(M == 1, ties.method = "first")
keep <- f < ncol(M)
M <- M[keep, , drop = FALSE]
f <- f[keep]
cov <- cov[keep, , drop = FALSE]

sex <- as.integer(cov$Sex)
sex[is.na(sex)] <- 3L # Unknown, matching the IPM
tag <- as.integer(cov$Tag_type)
nind <- nrow(M)
nocc <- ncol(M)

cat(sprintf(
  "CJS data: %d individuals, %d occasions (%d-%d)\n",
  nind, nocc, min(YEARS), max(YEARS)
))
cat("tag:", paste(c("untagged", "geo", "GPS"), table(tag), sep = "=", collapse = "  "), "\n")
cat("sex:", paste(c("F", "M", "Unk"), table(sex), sep = "=", collapse = "  "), "\n\n")

## ===========================================================================
## SECTION 2 -- MARGINALISED CJS LIKELIHOOD
## ===========================================================================
# Individual-specific phi and p, latent states integrated out analytically.
# Validated in scripts/LEYE_juv90_module.R against brute-force enumeration.
dCJSi <- nimbleFunction(
  run = function(
    x = double(1), phi = double(1), p = double(1),
    f = double(0), last = double(0), log = integer(0, default = 0)
  ) {
    returnType(double(0))
    L <- 0
    for (t in f:last) {
      if (x[t] == 1) L <- t
    }
    # chi[t] = P(no detection in (t, last] | alive at t)
    chiF <- nimNumeric(1, value = 1, length = 1)
    ch <- 1
    if (last > f) {
      for (k in 1:(last - f)) {
        t <- last - k
        ch <- (1 - phi[t]) + phi[t] * (1 - p[t + 1]) * ch
      }
    }
    chi_f <- ch
    if (L <= f) {
      # released at f, never seen again
      if (log) return(log(chi_f))
      return(chi_f)
    }
    val <- 1
    for (t in (f + 1):L) {
      if (x[t] == 1) val <- val * p[t] else val <- val * (1 - p[t])
    }
    for (t in f:(L - 1)) val <- val * phi[t]
    chL <- 1
    if (L < last) {
      for (k in 1:(last - L)) {
        t <- last - k
        chL <- (1 - phi[t]) + phi[t] * (1 - p[t + 1]) * chL
      }
    }
    out <- val * chL
    if (log) return(log(out))
    return(out)
  }
)
rCJSi <- nimbleFunction(
  run = function(
    n = integer(0), phi = double(1), p = double(1),
    f = double(0), last = double(0)
  ) {
    returnType(double(1))
    return(nimNumeric(length = 17))
  }
)
registerDistributions(list(dCJSi = list(
  BUGSdist = "dCJSi(phi, p, f, last)",
  types = c("value = double(1)", "phi = double(1)", "p = double(1)"),
  discrete = TRUE
)))

## ===========================================================================
## SECTION 3 -- MODELS
## ===========================================================================
# durMat[i,t] = years since banding at the START of interval t (0 for untagged),
# centred so the intercept keeps its "effect at deployment" meaning.
durMat <- matrix(0, nind, nocc - 1)
for (i in 1:nind) {
  for (t in f[i]:(nocc - 1)) durMat[i, t] <- t - f[i]
}
durMat[tag == 1, ] <- 0 # slope applies to tagged birds only

# firstMat[i,t] = 1 if t is the first interval after banding
firstMat <- matrix(0, nind, nocc - 1)
for (i in 1:nind) firstMat[i, f[i]] <- 1

make_code <- function(transience = FALSE, duration = FALSE) {
  nimbleCode({
    mu_phi ~ dnorm(0, sd = 1.5)
    sd_phi ~ dexp(2)
    for (t in 1:(nocc - 1)) {
      eps[t] ~ dnorm(0, sd = sd_phi)
      lphi[t] <- mu_phi + eps[t]
    }
    bSex[1] <- 0
    bTag[1] <- 0
    for (s in 2:3) {
      bSex[s] ~ dnorm(0, sd = 1.5)
    }
    for (g in 2:3) {
      bTag[g] ~ dnorm(0, sd = 1.5)
    }
    lp ~ dnorm(0, sd = 1.5)
    bTagP[1] <- 0
    for (g in 2:3) {
      bTagP[g] ~ dnorm(0, sd = 1.5)
    }
    bTrans ~ dnorm(0, sd = 1.5) # first-interval offset
    bDur ~ dnorm(0, sd = 1) # per-year change in tag effect

    for (i in 1:nind) {
      logit(pind[i]) <- lp + bTagP[tag[i]]
      for (t in 1:nocc) {
        pvec[i, t] <- pind[i]
      }
      for (t in 1:(nocc - 1)) {
        logit(phiM[i, t]) <- lphi[t] + bSex[sex[i]] + bTag[tag[i]] +
          TRANS * bTrans * firstM[i, t] +
          DUR * bDur * durM[i, t]
      }
      y[i, 1:nocc] ~ dCJSi(phiM[i, 1:(nocc - 1)], pvec[i, 1:nocc], fcap[i], nocc)
    }
  })
}

fit_one <- function(label, TRANS, DUR) {
  code <- make_code()
  consts <- list(
    nind = nind, nocc = nocc, sex = sex, tag = tag, fcap = f,
    firstM = firstMat, durM = durMat, TRANS = TRANS, DUR = DUR
  )
  inits <- function() {
    list(
      mu_phi = qlogis(.65), sd_phi = .3, eps = rnorm(nocc - 1, 0, .2),
      bSex = c(NA, 0, 0), bTag = c(NA, 0, 0), lp = 0, bTagP = c(NA, 0, 0),
      bTrans = 0, bDur = 0
    )
  }
  mons <- c("mu_phi", "sd_phi", "bSex", "bTag", "lp", "bTagP")
  if (TRANS) mons <- c(mons, "bTrans")
  if (DUR) mons <- c(mons, "bDur")
  Rm <- nimbleModel(code, constants = consts, data = list(y = M), inits = inits())
  conf <- configureMCMC(Rm, monitors = mons, print = FALSE, enableWAIC = TRUE)
  mc <- buildMCMC(conf)
  Cm <- compileNimble(Rm)
  Cmc <- compileNimble(mc, project = Rm)
  out <- runMCMC(Cmc,
    niter = 30000, nburnin = 10000, thin = 5, nchains = 3,
    inits = inits, samplesAsCodaMCMC = TRUE, setSeed = 1:3, WAIC = TRUE
  )
  cat(sprintf("\n--- %s | WAIC = %.1f ---\n", label, out$WAIC$WAIC))
  list(samples = out$samples, waic = out$WAIC$WAIC, label = label)
}

if (!exists("SKIP_TAG_RUN")) SKIP_TAG_RUN <- FALSE
if (!SKIP_TAG_RUN) {
  m0 <- fit_one("M0  sex + tag", FALSE, FALSE)
  m1 <- fit_one("M1  + transience", TRUE, FALSE)
  m2 <- fit_one("M2  + transience + tag x duration", TRUE, TRUE)

  orq <- function(s, cn) {
    x <- as.matrix(s)[, cn]
    sprintf(
      "%.2f (%.2f-%.2f) P(<1)=%.3f", exp(mean(x)),
      exp(quantile(x, .025)), exp(quantile(x, .975)), mean(x < 0)
    )
  }
  cat("\n=====================================================================\n")
  cat("TAG EFFECTS ON SURVIVAL (odds ratios vs untagged)\n")
  cat("=====================================================================\n")
  for (m in list(m0, m1, m2)) {
    cat(sprintf(
      "%-34s  geo %s | GPS %s\n", m$label,
      orq(m$samples, "bTag[2]"), orq(m$samples, "bTag[3]")
    ))
  }
  cat("\nWAIC:", sprintf("%s=%.1f", c("M0", "M1", "M2"), c(m0$waic, m1$waic, m2$waic)), "\n")

  q <- function(s, cn) {
    x <- as.matrix(s)[, cn]
    sprintf("%.3f (%.3f, %.3f)", mean(x), quantile(x, .025), quantile(x, .975))
  }
  cat("\ntransience (logit offset, first interval): ", q(m1$samples, "bTrans"), "\n")
  cat("duration slope (per yr since banding)    : ", q(m2$samples, "bDur"), "\n")
  cat("  P(duration slope < 0):", round(mean(as.matrix(m2$samples)[, "bDur"] < 0), 3), "\n")

  saveRDS(list(m0 = m0, m1 = m1, m2 = m2), "LEYE_tag_effects_CJS.rds")
  cat("\nWrote: LEYE_tag_effects_CJS.rds\n")
}
