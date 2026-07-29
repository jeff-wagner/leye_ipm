# =============================================================================
# Controlled comparison: LINEAR year vs FREE (factor) index
# All 3 weeks of data (n = 855). Habitat covariates EXCLUDED.
#
# Both models are identical except for the year term:
#   Model L (linear): log(mu) = a0 + a1*year_c + weekFE + b_et*log_et
#                               + b_ed*log_ed + eps_wet     [single slope a1]
#   Model F (factor): log(mu) = b0 + yearFE[t] + weekFE + b_et*log_et
#                               + b_ed*log_ed + eps_wet     [12 free year effects]
#
# Shared: negative binomial (var = mu + mu^2/r), week as nuisance factor,
#         centered log1p effort, non-centered wetland random intercept.
# This isolates the cost of the straight-line trend assumption.
# =============================================================================

library(nimble)
library(coda)
set.seed(20260521)

# ---- 1. Data prep (shared by both models) ----------------------------------
d <- read.csv("./Arin LEYE R code/data_2025.csv")
d$yr_idx <- as.integer(factor(d$year)) # 1..12
d$wk_idx <- as.integer(factor(d$week)) # 1..3
d$wet_idx <- as.integer(factor(d$wetland)) # 1..32

n <- nrow(d)
nYear <- length(unique(d$yr_idx))
nWeek <- length(unique(d$wk_idx))
nWet <- length(unique(d$wet_idx))

log_et <- log1p(d$et)
log_et <- log_et - mean(log_et)
log_ed <- log1p(d$ed)
log_ed <- log_ed - mean(log_ed)

# Centered calendar year (1 unit = 1 yr), both per-observation and per-level
yr_levels <- sort(unique(d$year))
yr_mean <- mean(yr_levels)
year_c <- d$year - yr_mean # per observation (Model L)
yr_c <- yr_levels - yr_mean # per year level (Model F slope)

# ---- 2. Model L: linear year -----------------------------------------------
codeL <- nimbleCode({
  a0 ~ dnorm(0, sd = 5)
  a1 ~ dnorm(0, sd = 5)
  b_et ~ dnorm(0, sd = 5)
  b_ed ~ dnorm(0, sd = 5)
  weekFE[1] <- 0
  for (w in 2:nWeek) {
    weekFE[w] ~ dnorm(0, sd = 5)
  }

  # --- ADD: free year deviations around the line (David's ey[t]) ---
  sd_ey ~ dexp(1)
  for (t in 1:nYear) {
    ey[t] ~ dnorm(0, sd = sd_ey)
  }
  # ----------------------------------------------------------------

  sd_wet ~ dexp(1)
  for (j in 1:nWet) {
    z_wet[j] ~ dnorm(0, sd = 1)
    eps_wet[j] <- z_wet[j] * sd_wet
  }
  r ~ dexp(0.1)
  for (i in 1:n) {
    log(mu[i]) <- a0 +
      a1 * year_c[i] +
      ey[yr[i]] +
      weekFE[wk[i]] +
      b_et * log_et[i] +
      b_ed * log_ed[i] +
      eps_wet[wet[i]]
    p[i] <- r / (r + mu[i])
    count[i] ~ dnbinom(prob = p[i], size = r)
  }
  for (t in 1:nYear) {
    lineL[t] <- exp(a0 + a1 * yr_c[t]) # pure linear trend
    levelL[t] <- exp(a0 + a1 * yr_c[t] + ey[t]) # trend + year deviation
    sd_ey_out <- sd_ey # convenience monitor
  }
  pctL <- 100 * (exp(a1) - 1) # % change / yr
})

constL <- list(
  n = n,
  nYear = nYear,
  nWeek = nWeek,
  nWet = nWet,
  yr = d$yr_idx,
  wk = d$wk_idx,
  wet = d$wet_idx,
  year_c = year_c,
  yr_c = yr_c,
  log_et = log_et,
  log_ed = log_ed
)
dataL <- list(count = d$count)
initsL <- function() {
  list(
    a0 = log(mean(d$count) + 0.1),
    a1 = 0,
    b_et = 0,
    b_ed = 0,
    weekFE = c(NA, rnorm(nWeek - 1, 0, .3)),
    ey = rnorm(nYear, 0, .3), # <-- ADD
    sd_ey = .3, # <-- ADD
    sd_wet = .5,
    z_wet = rnorm(nWet, 0, 1),
    r = 1
  )
}
monL <- c(
  "a0",
  "a1",
  "b_et",
  "b_ed",
  "weekFE",
  "sd_ey",
  "sd_wet",
  "r",
  "lineL",
  "levelL",
  "pctL"
)

# ---- 3. Model F: free factor index -----------------------------------------
codeF <- nimbleCode({
  b0 ~ dnorm(0, sd = 5)
  b_et ~ dnorm(0, sd = 5)
  b_ed ~ dnorm(0, sd = 5)
  yearFE[1] <- 0
  for (t in 2:nYear) {
    yearFE[t] ~ dnorm(0, sd = 5)
  }
  weekFE[1] <- 0
  for (w in 2:nWeek) {
    weekFE[w] ~ dnorm(0, sd = 5)
  }
  sd_wet ~ dexp(1)
  for (j in 1:nWet) {
    z_wet[j] ~ dnorm(0, sd = 1)
    eps_wet[j] <- z_wet[j] * sd_wet
  }
  r ~ dexp(0.1)
  for (i in 1:n) {
    log(mu[i]) <- b0 +
      yearFE[yr[i]] +
      weekFE[wk[i]] +
      b_et * log_et[i] +
      b_ed * log_ed[i] +
      eps_wet[wet[i]]
    p[i] <- r / (r + mu[i])
    count[i] ~ dnbinom(prob = p[i], size = r)
  }
  for (t in 1:nYear) {
    index[t] <- exp(b0 + yearFE[t])
  } # free index, ref week
  meanFE <- mean(yearFE[1:nYear])
  for (t in 1:nYear) {
    num_t[t] <- yr_c[t] * (yearFE[t] - meanFE)
    den_t[t] <- yr_c[t] * yr_c[t]
  }
  b_trend <- sum(num_t[1:nYear]) / sum(den_t[1:nYear])
  pctF <- 100 * (exp(b_trend) - 1) # % change / yr
})

constF <- list(
  n = n,
  nYear = nYear,
  nWeek = nWeek,
  nWet = nWet,
  yr = d$yr_idx,
  wk = d$wk_idx,
  wet = d$wet_idx,
  yr_c = yr_c,
  log_et = log_et,
  log_ed = log_ed
)
dataF <- list(count = d$count)
initsF <- function() {
  list(
    b0 = log(mean(d$count) + 0.1),
    b_et = 0,
    b_ed = 0,
    yearFE = c(NA, rnorm(nYear - 1, 0, .3)),
    weekFE = c(NA, rnorm(nWeek - 1, 0, .3)),
    sd_wet = .5,
    z_wet = rnorm(nWet, 0, 1),
    r = 1
  )
}
monF <- c(
  "b0",
  "b_et",
  "b_ed",
  "yearFE",
  "weekFE",
  "sd_wet",
  "r",
  "index",
  "b_trend",
  "pctF"
)

# ---- 4. Run -----------------------------------------------------------------
run_nimble <- function(
  code,
  const,
  data,
  inits,
  mon,
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
  Rc <- buildMCMC(configureMCMC(Rm, monitors = mon))
  Cm <- compileNimble(Rm)
  Cc <- compileNimble(Rc, project = Rm)
  runMCMC(
    Cc,
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    nchains = nchains,
    inits = inits,
    samplesAsCodaMCMC = TRUE,
    setSeed = 1:nchains
  )
}

cat("=== Model L (linear year) ===\n")
sampL <- run_nimble(codeL, constL, dataL, initsL, monL)
cat("=== Model F (factor index) ===\n")
sampF <- run_nimble(codeF, constF, dataF, initsF, monF)
matL <- as.matrix(sampL)
matF <- as.matrix(sampF)

# ---- 5. Convergence ---------------------------------------------------------
cat("\n--- Rhat (want < 1.05) ---\n")
print(round(
  gelman.diag(sampL[, c("a0", "a1", "sd_wet", "r")], multivariate = FALSE)$psrf,
  3
))
print(round(
  gelman.diag(
    sampF[, c("b0", "b_trend", "pctF", "sd_wet", "r")],
    multivariate = FALSE
  )$psrf,
  3
))

# ---- 6. Slopes side by side -------------------------------------------------
q <- function(x) c(mean(x), quantile(x, c(.025, .975)))
cat("\n=========== OVERALL TREND (% change / yr) ===========\n")
cat(sprintf(
  "Model L (linear slope):    %6.2f%%  (95%% CrI %6.2f, %6.2f)  P(decline)=%.3f\n",
  q(matL[, "pctL"])[1],
  q(matL[, "pctL"])[2],
  q(matL[, "pctL"])[3],
  mean(matL[, "a1"] < 0)
))
cat(sprintf(
  "Model F (factor-derived):  %6.2f%%  (95%% CrI %6.2f, %6.2f)  P(decline)=%.3f\n",
  q(matF[, "pctF"])[1],
  q(matF[, "pctF"])[2],
  q(matF[, "pctF"])[3],
  mean(matF[, "b_trend"] < 0)
))
cat("=====================================================\n")

# ---- 7. Relative-scale trajectories for overlay ----------------------------
rel_summ <- function(mat, cols) {
  X <- mat[, cols, drop = FALSE]
  Xr <- X / rowMeans(X)
  data.frame(
    year = yr_levels,
    mean = colMeans(Xr),
    lwr = apply(Xr, 2, quantile, .025),
    upr = apply(Xr, 2, quantile, .975)
  )
}
F_index <- rel_summ(matF, paste0("index[", 1:nYear, "]"))
L_line <- rel_summ(matL, paste0("lineL[", 1:nYear, "]")) # pure linear trend
L_level <- rel_summ(matL, paste0("levelL[", 1:nYear, "]")) # trend + year deviation

# ---- 8. Overlay plot --------------------------------------------------------
png("leye_linear_vs_factor_allData.png", width = 1700, height = 1100, res = 200)
par(mar = c(4.5, 4.5, 2.5, 1))
yr <- yr_levels
ylim <- range(c(F_index$lwr, F_index$upr, L_line$lwr, L_line$upr, L_level$mean))
plot(
  yr,
  F_index$mean,
  type = "n",
  ylim = ylim,
  xlab = "Year",
  ylab = "Relative abundance index (mean year = 1)",
  main = "All 3 weeks: free annual index vs linear trend (NB, habitat excluded)"
)
abline(h = 1, col = "grey80", lty = 3)
polygon(
  c(yr, rev(yr)),
  c(L_line$lwr, rev(L_line$upr)),
  col = adjustcolor("darkorange", 0.15),
  border = NA
)
lines(yr, L_line$mean, col = "darkorange", lwd = 2.5)
# Model L realized yearly level (trend + ey): open circles, shows what ey absorbs
points(yr, L_level$mean, pch = 1, col = "darkorange", cex = 1.1)
arrows(
  yr,
  F_index$lwr,
  yr,
  F_index$upr,
  length = .03,
  angle = 90,
  code = 3,
  col = "steelblue"
)
points(yr, F_index$mean, pch = 16, col = "steelblue", cex = 1.2)
lines(yr, F_index$mean, col = "steelblue", lwd = 1, lty = 2)
legend(
  "topright",
  bty = "n",
  cex = .85,
  legend = c(
    "Model F: free annual index +/- 95% CrI",
    "Model L: pure linear trend +/- 95% CrI",
    "Model L: realized level (incl. year RE)"
  ),
  col = c("steelblue", "darkorange", "darkorange"),
  pch = c(16, NA, 1),
  lty = c(2, 1, NA),
  lwd = c(1, 2.5, NA)
)
mtext(
  sprintf(
    "Slope L = %.1f%%/yr | Slope F = %.1f%%/yr",
    q(matL[, "pctL"])[1],
    q(matF[, "pctF"])[1]
  ),
  side = 1,
  line = -1.1,
  cex = .8
)
dev.off()
cat("\nWrote: leye_linear_vs_factor_3wk.png\n")

# Overlay plot Week 3 (Arin/David model) vs. All Data model
png(
  "leye_linear_wk3_vs_factor_allData.png",
  width = 1700,
  height = 1100,
  res = 200
)
par(mar = c(4.5, 4.5, 2.5, 1))
yr <- yr_levels
ylim <- range(c(F_index$lwr, F_index$upr, L_line$lwr, L_line$upr, L_level$mean))
plot(
  yr,
  F_index$mean,
  type = "n",
  ylim = ylim,
  xlab = "Year",
  ylab = "Relative abundance index (mean year = 1)",
  main = "Linear trend (week 3 only, habitat included) vs\nFree annual index (all data, habitat excluded)"
)
abline(h = 1, col = "grey80", lty = 3)
polygon(
  c(yr, rev(yr)),
  c(A_line$lwr, rev(A_line$upr)),
  col = adjustcolor("forestgreen", .15),
  border = NA
)
lines(yr, A_line$mean, col = "forestgreen", lwd = 2.5) # matched linear control
arrows(
  yr,
  F_index$lwr,
  yr,
  F_index$upr,
  length = .03,
  angle = 90,
  code = 3,
  col = "steelblue"
)
points(yr, F_index$mean, pch = 16, col = "steelblue", cex = 1.2)
lines(yr, F_index$mean, col = "steelblue", lwd = 1, lty = 2)
legend(
  "topright",
  bty = "n",
  cex = .85,
  legend = c(
    "Model A: linear trend (incl. year RE) +/- 95% CrI",
    "Model B: free annual index +/- 95% CrI"
  ),
  col = c("forestgreen", "steelblue"),
  pch = c(NA, 16),
  lty = c(1, 2),
  lwd = c(2.5, 1)
)
mtext(
  sprintf(
    "Slope A = %.1f%%/yr | Slope B = %.1f%%/yr",
    q(A_pct)[1],
    q(matF[, "pctF"])[1]
  ),
  side = 1,
  line = -1.1,
  cex = .8
)
dev.off()
cat("\nWrote: leye_linear_wk3_vs_factor_allData.png\n")

# ---- 9. Curvature check -----------------------------------------------------
dev_from_line <- F_index$mean - L_line$mean
cat("\nPer-year deviation of free index from linear trend (relative units):\n")
print(round(setNames(dev_from_line, yr), 3))
cat(sprintf(
  "RMS deviation: %.3f  (large => straight line is smoothing over real structure)\n",
  sqrt(mean(dev_from_line^2))
))
