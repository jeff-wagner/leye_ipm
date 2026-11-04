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
# Transience: untagged birds carry about twice the share of first-intervals
# (16.6%) that geolocator birds do (8.4%; GPS 13.8%), because tagged cohorts
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
#
# ALSO REPORTED (Section 1B, no model fitting): the exposure imbalance that
# motivates M1, sex composition and banding-era overlap by tag class, and the
# mass-at-banding comparison. Those numbers appear in the manuscript as
# confounds to disclose, and were previously computed at the console with no
# record -- so when the encounter history changed they could not be re-derived.
# They are now printed from the same objects the models are fitted to.
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
## SECTION 1B -- EXPOSURE AND COVARIATE CONFOUNDS (no model fitting)
## ===========================================================================
# Every tag-related number the manuscript quotes outside the MCMC used to be
# computed by hand at the console, which is why the morphometric comparison had
# no script and could not be re-derived when the encounter history changed.
# These are printed from the same objects the models are fitted to, so the
# analysis set is identical by construction.
TAGLAB <- c("untagged", "geolocator", "GPS")

# Exposure: share of each tag class's intervals that are FIRST intervals after
# banding. This is the imbalance that motivates M1 -- if survival differs in the
# first interval, the class carrying more of them is affected most.
#
# TWO denominators, because the answer depends on which one and the manuscript
# needs the defensible one. "Monitored" counts every interval from banding to
# the end of the study: that is what the marginalised likelihood actually uses,
# since intervals after the last detection still contribute through chi. "To
# last detection" counts only intervals up to the final resighting. The second
# is the more intuitive reading of "exposure" and it REVERSES the ordering --
# GPS birds, being resighted least, have the largest first-interval share of it.
# The transience argument depends on the first: under the likelihood-relevant
# denominator the untagged class carries about twice the first-interval share of
# the geolocator class, which is the imbalance M1 controls for.
cat("--- exposure (first-interval share, two denominators) ---\n")
last_seen <- apply(M == 1, 1, function(r) max(which(r)))
for (g in 1:3) {
  i <- which(tag == g)
  monit <- sum(nocc - f[i]) # f .. nocc-1, every interval in the likelihood
  tolast <- sum(pmax(last_seen[i] - f[i], 0)) + length(i) # f .. last, min 1
  cat(sprintf(
    "  %-11s n=%3d  first=%3d | monitored=%4d (%.1f%%)  to-last-detection=%4d (%.1f%%)\n",
    TAGLAB[g], length(i), length(i), monit, 100 * length(i) / monit,
    tolast, 100 * length(i) / tolast
  ))
}

# Sex composition. Males survive better, so a female-skewed tag class is biased
# DOWNWARD; the sex covariate is in the model for this reason.
cat("\n--- sex composition by tag class (% of class) ---\n")
sx <- round(100 * prop.table(table(TAGLAB[tag], c("F", "M", "Unk")[sex]), 1), 1)
print(sx)

# Banding-era overlap. Tag type is largely confounded with period; what keeps the
# contrast from being purely across eras is the untagged birds banded in the same
# window as the GPS birds.
yb <- cov$yr_banded
cat("\n--- banding years by tag class ---\n")
for (g in 1:3) {
  i <- which(tag == g)
  cat(sprintf(
    "  %-11s %d-%d (median %d)\n", TAGLAB[g],
    min(yb[i], na.rm = TRUE), max(yb[i], na.rm = TRUE),
    median(yb[i], na.rm = TRUE)
  ))
}
gw <- range(yb[tag == 3], na.rm = TRUE)
cat(sprintf(
  "  untagged banded inside the GPS window (%d-%d): %d birds\n",
  gw[1], gw[2], sum(tag == 1 & yb >= gw[1] & yb <= gw[2], na.rm = TRUE)
))

# Selection bias on body mass at banding. If tagged birds were a heavier --
# i.e. better-conditioned -- subset, the tag effect is CONSERVATIVE, because
# condition should push their survival the other way.
cat("\n--- mass at banding (g), untagged as reference ---\n")
md <- data.frame(mass = cov$init_mass, grp = factor(TAGLAB[tag], levels = TAGLAB))
md <- md[is.finite(md$mass), ]
for (g in TAGLAB) {
  x <- md$mass[md$grp == g]
  cat(sprintf("  %-11s n=%3d  mean=%.1f  sd=%.1f\n", g, length(x), mean(x), sd(x)))
}
fitm <- lm(mass ~ grp, data = md)
ci <- confint(fitm)
sm <- summary(fitm)$coefficients
for (g in TAGLAB[-1]) {
  nm <- paste0("grp", g)
  cat(sprintf(
    "  %-11s vs untagged: %+.2f g (95%% CI %.2f to %.2f), p = %s\n",
    g, sm[nm, 1], ci[nm, 1], ci[nm, 2], format.pval(sm[nm, 4], digits = 3, eps = 1e-4)
  ))
}
cat("\n")

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

## ===========================================================================
## SECTION 4 -- FIT (or reload) AND REPORT
## ===========================================================================
# The reporting runs off the cache as readily as off a fresh fit:
#     SKIP_TAG_RUN <- TRUE; source("scripts/LEYE_tag_effects_CJS.R")
# reprints every number below from LEYE_tag_effects_CJS.rds in seconds. The
# three models take about twenty minutes to refit, and a quoted value that can
# only be re-derived by paying that is a value nobody re-derives -- which is
# how the previous set came to be quoted for three months against an encounter
# history that had since been corrected.
OUT <- "LEYE_tag_effects_CJS.rds"
if (!exists("SKIP_TAG_RUN")) SKIP_TAG_RUN <- FALSE

if (SKIP_TAG_RUN) {
  if (!file.exists(OUT)) stop("SKIP_TAG_RUN is set but ", OUT, " does not exist.")
  fits <- readRDS(OUT)
  m0 <- fits$m0
  m1 <- fits$m1
  m2 <- fits$m2
  cat(sprintf("Reporting from cached %s (fitted %s)\n", OUT,
              format(file.mtime(OUT), "%Y-%m-%d")))
} else {
  m0 <- fit_one("M0  sex + tag", FALSE, FALSE)
  m1 <- fit_one("M1  + transience", TRUE, FALSE)
  m2 <- fit_one("M2  + transience + tag x duration", TRUE, TRUE)
  saveRDS(list(m0 = m0, m1 = m1, m2 = m2), OUT)
  cat("\nWrote: ", OUT, "\n")
}

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
cat("  P(transience offset > 0):", round(mean(as.matrix(m1$samples)[, "bTrans"] > 0), 3), "\n")
cat("duration slope (per yr since banding)    : ", q(m2$samples, "bDur"), "\n")
cat("  P(duration slope < 0):", round(mean(as.matrix(m2$samples)[, "bDur"] < 0), 3), "\n")

# Convergence, on the parameters this script actually quotes. Reported rather
# than assumed: these are three-chain fits of 20000 post-burnin draws apiece,
# and the main IPM has already produced one companion whose Rhat sat at 1.19
# while every number drawn from it was tabulated as fact.
cat("\n--- Rhat on the quantities reported above ---\n")
for (m in list(m0, m1, m2)) {
  p <- intersect(c("bTag[2]", "bTag[3]", "bSex[2]", "bSex[3]", "bTrans", "bDur",
                   "mu_phi", "sd_phi", "lp"), colnames(m$samples[[1]]))
  r <- gelman.diag(m$samples[, p], multivariate = FALSE)$psrf[, 1]
  cat(sprintf("  %-34s max Rhat %.3f (%s)%s\n", m$label, max(r),
              names(which.max(r)), if (max(r) > 1.05) "   <- CHECK" else ""))
}
