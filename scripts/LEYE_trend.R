# =============================================================================
# Lesser Yellowlegs (LEYE) annual abundance index + overall trend
# Negative-binomial GLMM fit in NIMBLE
#
# Response : weekly survey counts (kept un-aggregated; 3 visits = replication)
# Index    : year as a factor  -> annual relative-abundance index w/ CIs
# Trend    : derived slope of year effects -> single "% change per year" w/ CI
# Nuisance : log effort (et, ed), week factor, wetland random intercept
#
# Model:  y ~ NB(mu, r),  var = mu + mu^2/r
#   log(mu) = b0 + yearFE[t] + weekFE[w] + b_et*log_et + b_ed*log_ed + eps_wet
# =============================================================================

library(nimble)
library(coda)

set.seed(20260521)

# ---- 1. Data prep -----------------------------------------------------------
d <- read.csv("./Arin LEYE R code/data_2025.csv")

# Integer indices (1..K) for factor-like terms
d$year_f <- factor(d$year)
d$yr_idx <- as.integer(d$year_f) # 1..12
d$wk_idx <- as.integer(factor(d$week)) # 1..3
d$wet_idx <- as.integer(factor(d$wetland)) # 1..32

nYear <- nlevels(d$year_f)
nWeek <- length(unique(d$wk_idx))
nWet <- length(unique(d$wet_idx))
N <- nrow(d)

# Centered log effort. log1p handles the zeros (ed==0 on ~29% of surveys; et has one 0).
# Centering -> b0/year effects interpret at mean effort, and improves mixing.
log_et <- log1p(d$et)
log_et <- log_et - mean(log_et)
log_ed <- log1p(d$ed)
log_ed <- log_ed - mean(log_ed)

# Numeric year, centered, for the derived trend (1 unit = 1 year)
yr_num_c <- as.numeric(levels(d$year_f)) # the 12 actual years
yr_num_c <- yr_num_c - mean(yr_num_c) # centered, length nYear

# ---- 2. NIMBLE model --------------------------------------------------------
code <- nimbleCode({
  # ----- Priors -----
  b0 ~ dnorm(0, sd = 5)
  b_et ~ dnorm(0, sd = 5)
  b_ed ~ dnorm(0, sd = 5)

  # Year factor: corner constraint (year 1 = 0) so b0 is the reference-year mean.
  yearFE[1] <- 0
  for (t in 2:nYear) {
    yearFE[t] ~ dnorm(0, sd = 5)
  }

  # Week factor: corner constraint (week 1 = 0)
  weekFE[1] <- 0
  for (w in 2:nWeek) {
    weekFE[w] ~ dnorm(0, sd = 5)
  }

  # Wetland random intercept
  sd_wet ~ dexp(1) # weakly-informative, > 0
  for (j in 1:nWet) {
    # non-centered parameterization of the wetland random effect;
    z_wet[j] ~ dnorm(0, sd = 1) # breaks the funnel-style correlation between b0 and the RE
    eps_wet[j] <- z_wet[j] * sd_wet
  }

  # NB dispersion (mean/var form): var = mu + mu^2 / r
  r ~ dexp(0.1) # prior mean 10; allows wide range

  # ----- Likelihood -----
  for (i in 1:N) {
    log(mu[i]) <- b0 +
      yearFE[yr[i]] +
      weekFE[wk[i]] +
      b_et * log_et[i] +
      b_ed * log_ed[i] +
      eps_wet[wet[i]]
    p[i] <- r / (r + mu[i]) # convert mean -> NB prob
    y[i] ~ dnbinom(prob = p[i], size = r)
  }

  # ----- Derived quantities -----
  # (a) Annual index on the response scale: expected count at reference week (1)
  #     and mean effort (covariates centered -> 0), averaged over wetland RE (=0 at mean).
  for (t in 1:nYear) {
    index[t] <- exp(b0 + yearFE[t])
  }

  # (b) Overall trend = OLS slope of yearFE on centered numeric year, computed
  #     each iteration so its posterior CI is fully propagated.
  #     slope is on the log scale -> exp(slope) = multiplicative change per year.
  ybar <- 0 # yearFE already (effectively) centered via corner constraint; recompute mean
  meanFE <- mean(yearFE[1:nYear])
  for (t in 1:nYear) {
    num_t[t] <- yr_num_c[t] * (yearFE[t] - meanFE)
    den_t[t] <- yr_num_c[t] * yr_num_c[t]
  }
  b_trend <- sum(num_t[1:nYear]) / sum(den_t[1:nYear]) # log-scale slope/yr
  pct_change_yr <- 100 * (exp(b_trend) - 1) # % change per year
})

# ---- 3. Constants, data, inits ---------------------------------------------
constants <- list(
  N = N,
  nYear = nYear,
  nWeek = nWeek,
  nWet = nWet,
  yr = d$yr_idx,
  wk = d$wk_idx,
  wet = d$wet_idx,
  log_et = log_et,
  log_ed = log_ed,
  yr_num_c = yr_num_c
)

data <- list(y = d$count)

inits <- function() {
  list(
    b0 = log(mean(d$count) + 0.1),
    b_et = 0,
    b_ed = 0,
    yearFE = c(NA, rnorm(nYear - 1, 0, 0.3)),
    weekFE = c(NA, rnorm(nWeek - 1, 0, 0.3)),
    sd_wet = 0.5,
    eps_wet = rnorm(nWet, 0, 0.3),
    r = 1
  )
}

# ---- 4. Build, compile, run -------------------------------------------------
monitors <- c(
  "b0",
  "b_et",
  "b_ed",
  "yearFE",
  "weekFE",
  "sd_wet",
  "r",
  "index",
  "b_trend",
  "pct_change_yr"
)

Rmodel <- nimbleModel(
  code,
  constants = constants,
  data = data,
  inits = inits(),
  calculate = FALSE
)
conf <- configureMCMC(Rmodel, monitors = monitors)
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

samples <- runMCMC(
  Cmcmc,
  niter = 30000,
  nburnin = 10000,
  thin = 5,
  nchains = 3,
  inits = inits,
  samplesAsCodaMCMC = TRUE,
  setSeed = 1:3
)

# ---- 5. Diagnostics ---------------------------------------------------------
key <- c("b0", "b_et", "b_ed", "sd_wet", "r", "b_trend", "pct_change_yr")
print(summary(samples[, key]))
cat("\n--- Convergence (Gelman-Rubin Rhat; want < 1.05) ---\n")
gd <- gelman.diag(samples, multivariate = FALSE)$psrf
print(round(gd[c(key, paste0("index[", 1:nYear, "]")), ], 3))
cat("\n--- Effective sample size (key params) ---\n")
print(round(effectiveSize(samples[, key])))

# ---- 6. Summaries: trend + annual index ------------------------------------
mat <- as.matrix(samples)

qfun <- function(x) {
  c(
    mean = mean(x),
    lwr = quantile(x, .025, names = FALSE),
    upr = quantile(x, .975, names = FALSE)
  )
}

# Overall trend
trend_pct <- qfun(mat[, "pct_change_yr"])
prob_decline <- mean(mat[, "b_trend"] < 0)
cat("\n=========================================================\n")
cat(sprintf(
  "OVERALL TREND: %.2f%% per year (95%% CI %.2f to %.2f)\n",
  trend_pct["mean"],
  trend_pct["lwr.2.5%"],
  trend_pct["upr.97.5%"]
))
cat(sprintf("Posterior prob. of decline (slope < 0): %.3f\n", prob_decline))
cat("=========================================================\n")

# Annual index table
idx_cols <- paste0("index[", 1:nYear, "]")
index_summ <- t(apply(mat[, idx_cols], 2, qfun))
index_df <- data.frame(
  year = as.numeric(levels(d$year_f)),
  mean = index_summ[, 1],
  lwr = index_summ[, 2],
  upr = index_summ[, 3]
)
cat(
  "\n--- Annual abundance index (expected count, ref week, mean effort) ---\n"
)
print(round(index_df, 3), row.names = FALSE)

# ---- 7. Plot: annual index with 95% CrI ------------------------------------
png("leye_annual_index.png", width = 1500, height = 1000, res = 200)
par(mar = c(4.5, 4.5, 2, 1))
with(index_df, {
  plot(
    year,
    mean,
    type = "n",
    ylim = range(c(lwr, upr)),
    xlab = "Year",
    ylab = "Relative abundance index (expected count)",
    main = "LEYE annual abundance index (NB-GLMM, NIMBLE)"
  )
  polygon(
    c(year, rev(year)),
    c(lwr, rev(upr)),
    col = adjustcolor("steelblue", 0.20),
    border = NA
  )
  lines(year, mean, col = "steelblue", lwd = 2)
  points(year, mean, pch = 16, col = "steelblue")
})
# overlay trend direction text
mtext(
  sprintf(
    "Overall trend: %.1f%%/yr (95%% CrI %.1f to %.1f); P(decline)=%.2f",
    trend_pct["mean"],
    trend_pct["lwr.2.5%"],
    trend_pct["upr.97.5%"],
    prob_decline
  ),
  side = 3,
  line = -1.2,
  cex = 0.8
)
dev.off()

cat("\nWrote: leye_annual_index.png\n")

# ---- 8. (Optional) posterior predictive check on zero proportion -----------
# Quick check that the NB captures the 59% zeros. Uncomment to run.
pzero_obs <- mean(d$count == 0)
nsim <- 500
idx <- sample(nrow(mat), nsim)
pzero_rep <- numeric(nsim)
for (s in seq_len(nsim)) {
  ps <- mat[idx[s], ]
  mu_s <- exp(
    ps["b0"] +
      ps[paste0("yearFE[", d$yr_idx, "]")] +
      ps[paste0("weekFE[", d$wk_idx, "]")] +
      ps["b_et"] * log_et +
      ps["b_ed"] * log_ed +
      rnorm(nWet, 0, ps["sd_wet"])[d$wet_idx]
  )
  rr <- ps["r"]
  pp <- rr / (rr + mu_s)
  yrep <- rnbinom(N, size = rr, prob = pp)
  pzero_rep[s] <- mean(yrep == 0)
}
cat(sprintf(
  "PPC zero-prop: obs=%.3f, rep mean=%.3f (Bayesian p=%.3f)\n",
  pzero_obs,
  mean(pzero_rep),
  mean(pzero_rep >= pzero_obs)
))
