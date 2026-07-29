# =============================================================================
# Balanced-panel sensitivity check for the LEYE decline
#
# Question: is the ~ -8%/yr trend an artifact of changing site composition as
# the survey program grew (13 -> ~30 wetlands)? Refit the free-index NB model
# (Model F) on only the CORE wetlands surveyed in nearly every year, so the
# trend cannot be driven by which sites entered over time.
#
# Compares the balanced-panel slope + annual index to the full-data result.
# Identical model structure to Model F in leye_linear_vs_factor_3wk.R
# (NB, year factor, week nuisance, log effort, non-centered wetland RE).
# =============================================================================

library(nimble)
library(coda)
set.seed(20260521)

# ---- 0. Threshold: keep wetlands surveyed in >= MIN_YEARS distinct years ---
#   12 = strict (9 wetlands, 306 surveys);  11 = recommended (12 wetlands, 396)
MIN_YEARS <- 11

# ---- 1. Data prep -----------------------------------------------------------
d_full <- read.csv("./Arin LEYE R code/data_2025.csv")

yrs_per_wet <- tapply(d_full$year, d_full$wetland, function(x) {
  length(unique(x))
})
core_wet <- as.integer(names(yrs_per_wet[yrs_per_wet >= MIN_YEARS]))
cat(sprintf(
  "Balanced panel: %d wetlands with >= %d yrs of data (of %d total).\n",
  length(core_wet),
  MIN_YEARS,
  length(yrs_per_wet)
))

d <- d_full[d_full$wetland %in% core_wet, ]
cat(sprintf(
  "Surveys retained: %d of %d (%.0f%%)\n",
  nrow(d),
  nrow(d_full),
  100 * nrow(d) / nrow(d_full)
))

d$yr_idx <- as.integer(factor(d$year)) # 1..12 (all years still present)
d$wk_idx <- as.integer(factor(d$week)) # 1..3
d$wet_idx <- as.integer(factor(d$wetland)) # 1..nWet (re-indexed contiguous)

n <- nrow(d)
nYear <- length(unique(d$yr_idx))
nWeek <- length(unique(d$wk_idx))
nWet <- length(unique(d$wet_idx))
stopifnot(nYear == 12) # need all years for the index

log_et <- log1p(d$et)
log_et <- log_et - mean(log_et)
log_ed <- log1p(d$ed)
log_ed <- log_ed - mean(log_ed)
yr_levels <- sort(unique(d$year))
yr_c <- yr_levels - mean(yr_levels)

# ---- 2. Model F (free factor index), same structure as before --------------
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
  }
  meanFE <- mean(yearFE[1:nYear])
  for (t in 1:nYear) {
    num_t[t] <- yr_c[t] * (yearFE[t] - meanFE)
    den_t[t] <- yr_c[t] * yr_c[t]
  }
  b_trend <- sum(num_t[1:nYear]) / sum(den_t[1:nYear])
  pctF <- 100 * (exp(b_trend) - 1)
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

# ---- 3. Run -----------------------------------------------------------------
Rm <- nimbleModel(
  codeF,
  constants = constF,
  data = dataF,
  inits = initsF(),
  calculate = FALSE
)
Rc <- buildMCMC(configureMCMC(Rm, monitors = monF))
Cm <- compileNimble(Rm)
Cc <- compileNimble(Rc, project = Rm)
samp <- runMCMC(
  Cc,
  niter = 60000,
  nburnin = 20000,
  thin = 10,
  nchains = 3,
  inits = initsF,
  samplesAsCodaMCMC = TRUE,
  setSeed = 1:3
)
mat <- as.matrix(samp)

# ---- 4. Convergence ---------------------------------------------------------
cat("\n--- Rhat (want < 1.05) ---\n")
print(round(
  gelman.diag(
    samp[, c("b0", "b_trend", "pctF", "sd_wet", "r")],
    multivariate = FALSE
  )$psrf,
  3
))
cat(
  "Min ESS (index + slope):",
  round(min(effectiveSize(samp[, c(paste0("index[", 1:nYear, "]"), "pctF")]))),
  "\n"
)

# ---- 5. Slope --------------------------------------------------------------
q <- function(x) c(mean(x), quantile(x, c(.025, .975)))
cat("\n=========== BALANCED-PANEL TREND ===========\n")
cat(sprintf(
  "Wetlands >= %d yrs (n=%d sites, %d surveys)\n",
  MIN_YEARS,
  nWet,
  n
))
cat(sprintf(
  "Slope: %6.2f%%/yr  (95%% CrI %6.2f, %6.2f)  P(decline)=%.3f\n",
  q(mat[, "pctF"])[1],
  q(mat[, "pctF"])[2],
  q(mat[, "pctF"])[3],
  mean(mat[, "b_trend"] < 0)
))
cat("Compare to FULL-DATA Model F: about -8.4%/yr.\n")
cat("============================================\n")

# ---- 6. Annual index table + plot ------------------------------------------
idx_cols <- paste0("index[", 1:nYear, "]")
isum <- t(apply(mat[, idx_cols], 2, function(x) {
  c(mean(x), quantile(x, c(.025, .975)))
}))
index_df <- data.frame(
  year = yr_levels,
  mean = isum[, 1],
  lwr = isum[, 2],
  upr = isum[, 3]
)
cat(
  "\n--- Balanced-panel annual index (expected count, ref week, mean effort) ---\n"
)
print(round(index_df, 3), row.names = FALSE)

png("leye_balanced_panel_index.png", width = 1600, height = 1050, res = 200)
par(mar = c(4.5, 4.5, 2.5, 1))
with(index_df, {
  plot(
    year,
    mean,
    type = "n",
    ylim = range(c(lwr, upr)),
    xlab = "Year",
    ylab = "Relative abundance index (expected count)",
    main = sprintf(
      "LEYE balanced panel (wetlands surveyed >= %d yrs, n=%d sites)",
      MIN_YEARS,
      nWet
    )
  )
  polygon(
    c(year, rev(year)),
    c(lwr, rev(upr)),
    col = adjustcolor("seagreen", .18),
    border = NA
  )
  lines(year, mean, col = "seagreen", lwd = 2)
  points(year, mean, pch = 16, col = "seagreen")
})
mtext(
  sprintf(
    "Slope %.1f%%/yr (95%% CrI %.1f, %.1f); P(decline)=%.2f  | full-data was ~ -8.4%%/yr",
    q(mat[, "pctF"])[1],
    q(mat[, "pctF"])[2],
    q(mat[, "pctF"])[3],
    mean(mat[, "b_trend"] < 0)
  ),
  side = 3,
  line = -1.1,
  cex = .72
)
dev.off()
cat("\nWrote: leye_balanced_panel_index.png\n")
cat(
  "\nInterpretation: if this slope and index track the full-data result, the\n"
)
cat("decline is NOT an artifact of program expansion / site composition.\n")
