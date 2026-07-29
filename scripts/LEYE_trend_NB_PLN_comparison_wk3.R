# =============================================================================
# Three-model comparison on the SAME week-3 LEYE subset, with a CONTROLLED
# isolation of the trend term (linear vs factor) plus WAIC model selection.
#
#   Model A ("David"): Poisson-lognormal, LINEAR year (b0+b1*t+ey[t]),
#                      habitat covariates, wetland RE.        [whole approach]
#   Model B ("factor"): NB, year as a FACTOR (free index), log effort,
#                      wetland RE.                            [whole approach]
#   Model C ("linear"): NB, LINEAR year (b0+b1*yr), log effort, wetland RE,
#                      NO habitat, NO year RE.                [matches B except trend]
#
# Contrasts:
#   C vs B  -> isolates LINEAR vs FACTOR (everything else identical) <-- the test
#   A vs C  -> the remaining bundle (Poisson-LN vs NB; habitat; year RE)
#
# WAIC (conditional) is compared for B vs C: same RE structure, so the
# comparison is valid and answers "does the flexible year term earn its keep?"
# =============================================================================

library(nimble)
library(coda)
set.seed(20260521)

# ---- 1. Week-3 data prep (shared across all three) -------------------------
d <- read.csv("./Arin LEYE R code/data_2025.csv")
d3 <- d[d$week == 3, ]
d3$yr_idx <- as.integer(factor(d3$year))
d3$wet_idx <- as.integer(factor(d3$wetland))
n <- nrow(d3)
ny <- length(unique(d3$yr_idx))
nw <- length(unique(d3$wet_idx))

wet_tab <- d3[!duplicated(d3$wet_idx), ]
wet_tab <- wet_tab[order(wet_tab$wet_idx), ]
stopifnot(all(wet_tab$wet_idx == 1:nw))
wc <- as.numeric(scale(wet_tab$wc))
wa <- as.numeric(scale(wet_tab$wa))
wt <- as.integer(factor(wet_tab$wt))
nwt <- length(unique(wt))
wh <- as.integer(factor(wet_tab$wh))
nwh <- length(unique(wh))
wu <- as.integer(factor(wet_tab$wu))
nwu <- length(unique(wu))

et_z <- as.numeric(scale(d3$et))
ed_z <- as.numeric(scale(d3$ed))
log_et <- log1p(d3$et)
log_et <- log_et - mean(log_et)
log_ed <- log1p(d3$ed)
log_ed <- log_ed - mean(log_ed)
yr_c <- sort(unique(d3$year))
yr_c <- yr_c - mean(yr_c) # centered calendar year


# ---- 2. Shared runner (with WAIC) ------------------------------------------
run_nimble <- function(
  code,
  const,
  data,
  inits,
  mon,
  waic = TRUE,
  niter = 60000,
  nburnin = 20000,
  thin = 10,
  nchains = 3
) {
  Rm <- nimbleModel(
    code,
    constants = const,
    data = data,
    inits = inits(),
    calculate = FALSE
  )
  cf <- configureMCMC(Rm, monitors = mon, enableWAIC = waic)
  Rc <- buildMCMC(cf)
  Cm <- compileNimble(Rm)
  Cc <- compileNimble(Rc, project = Rm)
  out <- runMCMC(
    Cc,
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    nchains = nchains,
    inits = inits,
    samplesAsCodaMCMC = TRUE,
    WAIC = waic,
    setSeed = 1:nchains
  )
  if (waic) {
    list(samples = out$samples, waic = out$WAIC)
  } else {
    list(samples = out, waic = NULL)
  }
}
get_waic <- function(obj) tryCatch(obj$waic$WAIC, error = function(e) NA_real_)
get_pwaic <- function(obj) {
  tryCatch(obj$waic$pWAIC, error = function(e) NA_real_)
}

# ---- 3. Model A: David's Poisson-LN, linear + ey, habitat ------------------
# (tau.ey corrected from original 1/(sd.e*sd.ey) to 1/(sd.ey^2).)
codeA <- nimbleCode({
  for (i in 1:n) {
    count[i] ~ dpois(mu[i])
    log(mu[i]) <- bw[w[i]] + by[y[i]] + bet * et[i] + bed * ed[i] + e[i]
    e[i] ~ dnorm(0, tau.e)
  }
  for (j in 1:nw) {
    bw[j] <- bwc *
      wc[j] +
      bwa * wa[j] +
      bwt[wt[j]] +
      bwh[wh[j]] +
      bwu[wu[j]] +
      ew[j]
    ew[j] ~ dnorm(0, tau.ew)
  }
  for (t in 1:ny) {
    by[t] <- b0 + b1 * t + ey[t]
    ey[t] ~ dnorm(0, tau.ey)
  }
  bet ~ dnorm(0, 0.001)
  bed ~ dnorm(0, 0.001)
  bwc ~ dnorm(0, 0.001)
  bwa ~ dnorm(0, 0.001)
  bwt[1] <- 0
  bwh[1] <- 0
  bwu[1] <- 0
  for (k in 2:nwt) {
    bwt[k] ~ dnorm(0, 0.001)
  }
  for (k in 2:nwh) {
    bwh[k] ~ dnorm(0, 0.001)
  }
  for (k in 2:nwu) {
    bwu[k] ~ dnorm(0, 0.001)
  }
  b0 ~ dnorm(0, 0.001)
  b1 ~ dnorm(0, 0.001)
  tau.e <- 1 / (sd.e * sd.e)
  sd.e ~ dunif(0, 10)
  tau.ew <- 1 / (sd.ew * sd.ew)
  sd.ew ~ dunif(0, 10)
  tau.ey <- 1 / (sd.ey * sd.ey)
  sd.ey ~ dunif(0, 10)
  for (t in 1:ny) {
    lineA[t] <- exp(b0 + b1 * t)
  }
})
constA <- list(
  n = n,
  nw = nw,
  ny = ny,
  w = d3$wet_idx,
  y = d3$yr_idx,
  et = et_z,
  ed = ed_z,
  wc = wc,
  wa = wa,
  wt = wt,
  wh = wh,
  wu = wu,
  nwt = nwt,
  nwh = nwh,
  nwu = nwu
)
initsA <- function() {
  list(
    b0 = log(mean(d3$count) + 0.1),
    b1 = 0,
    bet = 0,
    bed = 0,
    bwc = 0,
    bwa = 0,
    bwt = c(NA, rep(0, nwt - 1)),
    bwh = c(NA, rep(0, nwh - 1)),
    bwu = c(NA, rep(0, nwu - 1)),
    ew = rnorm(nw, 0, .3),
    ey = rnorm(ny, 0, .3),
    e = rnorm(n, 0, .3),
    sd.e = .5,
    sd.ew = .5,
    sd.ey = .3
  )
}
monA <- c("b0", "b1", "bet", "bed", "sd.e", "sd.ew", "sd.ey", "lineA")

# ---- 4. Model B: NB factor (free index) ------------------------------------
codeB <- nimbleCode({
  b0 ~ dnorm(0, sd = 5)
  b_et ~ dnorm(0, sd = 5)
  b_ed ~ dnorm(0, sd = 5)
  yearFE[1] <- 0
  for (t in 2:ny) {
    yearFE[t] ~ dnorm(0, sd = 5)
  }
  sd_wet ~ dexp(1)
  for (j in 1:nw) {
    z_wet[j] ~ dnorm(0, sd = 1)
    eps_wet[j] <- z_wet[j] * sd_wet
  }
  r ~ dexp(0.1)
  for (i in 1:n) {
    log(mu[i]) <- b0 +
      yearFE[yr[i]] +
      b_et * log_et[i] +
      b_ed * log_ed[i] +
      eps_wet[wet[i]]
    p[i] <- r / (r + mu[i])
    count[i] ~ dnbinom(prob = p[i], size = r)
  }
  for (t in 1:ny) {
    index[t] <- exp(b0 + yearFE[t])
  }
  meanFE <- mean(yearFE[1:ny])
  for (t in 1:ny) {
    num_t[t] <- yr_c[t] * (yearFE[t] - meanFE)
    den_t[t] <- yr_c[t] * yr_c[t]
  }
  b_trend <- sum(num_t[1:ny]) / sum(den_t[1:ny])
  pct_change_yr <- 100 * (exp(b_trend) - 1)
})
constB <- list(
  n = n,
  ny = ny,
  nw = nw,
  yr = d3$yr_idx,
  wet = d3$wet_idx,
  log_et = log_et,
  log_ed = log_ed,
  yr_c = yr_c
)
initsB <- function() {
  list(
    b0 = log(mean(d3$count) + 0.1),
    b_et = 0,
    b_ed = 0,
    yearFE = c(NA, rnorm(ny - 1, 0, .3)),
    sd_wet = .5,
    z_wet = rnorm(nw, 0, 1),
    r = 1
  )
}
monB <- c(
  "b0",
  "b_et",
  "b_ed",
  "yearFE",
  "sd_wet",
  "r",
  "index",
  "b_trend",
  "pct_change_yr"
)

# ---- 5. Model C: NB linear -- matches B EXCEPT the year term ---------------
codeC <- nimbleCode({
  b0 ~ dnorm(0, sd = 5)
  b1 ~ dnorm(0, sd = 5)
  b_et ~ dnorm(0, sd = 5)
  b_ed ~ dnorm(0, sd = 5)
  sd_wet ~ dexp(1)
  for (j in 1:nw) {
    z_wet[j] ~ dnorm(0, sd = 1)
    eps_wet[j] <- z_wet[j] * sd_wet
  }
  r ~ dexp(0.1)
  for (i in 1:n) {
    log(mu[i]) <- b0 +
      b1 * yr_c[yr[i]] +
      b_et * log_et[i] +
      b_ed * log_ed[i] +
      eps_wet[wet[i]]
    p[i] <- r / (r + mu[i])
    count[i] ~ dnbinom(prob = p[i], size = r)
  }
  pct_change_yr <- 100 * (exp(b1) - 1)
  for (t in 1:ny) {
    lineC[t] <- exp(b0 + b1 * yr_c[t])
  }
})
constC <- constB # identical data inputs to B
initsC <- function() {
  list(
    b0 = log(mean(d3$count) + 0.1),
    b1 = 0,
    b_et = 0,
    b_ed = 0,
    sd_wet = .5,
    z_wet = rnorm(nw, 0, 1),
    r = 1
  )
}
monC <- c("b0", "b1", "b_et", "b_ed", "sd_wet", "r", "pct_change_yr", "lineC")

# ---- 6. Fit all three ------------------------------------------------------
cat("=== Model A (Poisson-LN, linear, habitat) ===\n")
fitA <- run_nimble(codeA, constA, list(count = d3$count), initsA, monA)
cat("=== Model B (NB, factor index) ===\n")
fitB <- run_nimble(codeB, constB, list(count = d3$count), initsB, monB)
cat("=== Model C (NB, linear; matches B except trend) ==\n")
fitC <- run_nimble(codeC, constC, list(count = d3$count), initsC, monC)
matA <- as.matrix(fitA$samples)
matB <- as.matrix(fitB$samples)
matC <- as.matrix(fitC$samples)

# ---- 7. Slopes side by side ------------------------------------------------
q <- function(x) c(mean(x), quantile(x, c(.025, .975)))
A_pct <- 100 * (exp(matA[, "b1"]) - 1)
C_pct <- matC[, "pct_change_yr"]
B_pct <- matB[, "pct_change_yr"]
cat("\n================= OVERALL TREND (% change / yr) =================\n")
cat(sprintf(
  "A  Poisson-LN linear : %6.2f%%  (%6.2f, %6.2f)  P(decl)=%.3f\n",
  q(A_pct)[1],
  q(A_pct)[2],
  q(A_pct)[3],
  mean(matA[, "b1"] < 0)
))
cat(sprintf(
  "C  NB linear (ctrl)  : %6.2f%%  (%6.2f, %6.2f)  P(decl)=%.3f\n",
  q(C_pct)[1],
  q(C_pct)[2],
  q(C_pct)[3],
  mean(matC[, "b1"] < 0)
))
cat(sprintf(
  "B  NB factor (deriv) : %6.2f%%  (%6.2f, %6.2f)  P(decl)=%.3f\n",
  q(B_pct)[1],
  q(B_pct)[2],
  q(B_pct)[3],
  mean(matB[, "b_trend"] < 0)
))
cat("-----------------------------------------------------------------\n")
cat(" C vs B  = pure linear-vs-factor effect (NB, no habitat, same RE)\n")
cat(" A vs C  = family (Pois-LN vs NB) + habitat + year-RE bundle\n")
cat("=================================================================\n")

# ---- 8. WAIC: does the factor model earn its extra parameters? -------------
wB <- get_waic(fitB)
wC <- get_waic(fitC)
wA <- get_waic(fitA)
cat("\n--- WAIC (lower = better predictive; conditional WAIC) ---\n")
cat(sprintf(
  "Model C (linear) WAIC = %8.2f   pWAIC = %6.2f\n",
  wC,
  get_pwaic(fitC)
))
cat(sprintf(
  "Model B (factor) WAIC = %8.2f   pWAIC = %6.2f\n",
  wB,
  get_pwaic(fitB)
))
cat(sprintf(
  "dWAIC (B - C) = %+.2f   [B & C matched except trend; valid comparison]\n",
  wB - wC
))
cat(sprintf(
  "Model A WAIC = %8.2f (NOT directly comparable: diff. family/RE/covars)\n",
  wA
))
dW <- wB - wC
verdict <- if (is.na(dW)) {
  "WAIC unavailable"
} else if (dW > 2) {
  "C (linear) preferred: the free year term does NOT earn its keep -> linearity is fine."
} else if (dW < -2) {
  "B (factor) preferred: the data want non-linear year structure -> linearity costs you."
} else {
  "Tie (|dWAIC|<2): no clear predictive gain from the factor model -> linearity defensible."
}
cat("VERDICT:", verdict, "\n")

# ---- 9. Visual control: B's free index vs C's matched linear line ----------
rel_summ <- function(mat, cols) {
  X <- mat[, cols, drop = FALSE]
  Xr <- X / rowMeans(X)
  data.frame(
    year = sort(unique(d3$year)),
    mean = colMeans(Xr),
    lwr = apply(Xr, 2, quantile, .025),
    upr = apply(Xr, 2, quantile, .975)
  )
}
B_index <- rel_summ(matB, paste0("index[", 1:ny, "]"))
C_line <- rel_summ(matC, paste0("lineC[", 1:ny, "]"))
A_line <- rel_summ(matA, paste0("lineA[", 1:ny, "]"))

png("leye_trend_isolation.png", width = 1700, height = 1100, res = 200)
par(mar = c(4.5, 4.5, 2.5, 1))
yr <- B_index$year
ylim <- range(c(B_index$lwr, B_index$upr, C_line$lwr, C_line$upr))
plot(
  yr,
  B_index$mean,
  type = "n",
  ylim = ylim,
  xlab = "Year",
  ylab = "Relative abundance index (mean year = 1)",
  main = "Linear vs factor year, isolated (both NB, same covariates & RE)"
)
abline(h = 1, col = "grey80", lty = 3)
polygon(
  c(yr, rev(yr)),
  c(C_line$lwr, rev(C_line$upr)),
  col = adjustcolor("forestgreen", .15),
  border = NA
)
lines(yr, C_line$mean, col = "forestgreen", lwd = 2.5) # matched linear control
lines(yr, A_line$mean, col = "darkorange", lwd = 1.5, lty = 2) # David's line, reference
arrows(
  yr,
  B_index$lwr,
  yr,
  B_index$upr,
  length = .03,
  angle = 90,
  code = 3,
  col = "steelblue"
)
points(yr, B_index$mean, pch = 16, col = "steelblue", cex = 1.2)
lines(yr, B_index$mean, col = "steelblue", lwd = 1, lty = 3)
legend(
  "topright",
  bty = "n",
  cex = .8,
  legend = c(
    "B: free annual index (NB) +/- 95% CrI",
    "C: matched linear trend (NB) +/- 95% CrI",
    "A: David linear (Pois-LN, habitat) ref"
  ),
  col = c("steelblue", "forestgreen", "darkorange"),
  pch = c(16, NA, NA),
  lty = c(3, 1, 2),
  lwd = c(1, 2.5, 1.5)
)
mtext(
  sprintf(
    "Slopes  A=%.1f  C=%.1f  B=%.1f  %%/yr   |   dWAIC(B-C)=%+.1f",
    q(A_pct)[1],
    q(C_pct)[1],
    q(B_pct)[1],
    dW
  ),
  side = 3,
  line = -1.1,
  cex = .75
)
dev.off()
cat("\nWrote: leye_trend_isolation.png\n")

# ---- 10. Numeric curvature: free index minus matched linear line -----------
dev_line <- B_index$mean - C_line$mean
cat(
  "\nPer-year deviation of free index from matched linear trend (relative units):\n"
)
print(round(setNames(dev_line, yr), 3))
cat(sprintf("RMS deviation: %.3f\n", sqrt(mean(dev_line^2))))
