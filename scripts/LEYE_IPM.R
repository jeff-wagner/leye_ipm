# =============================================================================
# Lesser Yellowlegs INTEGRATED POPULATION MODEL (IPM) -- NIMBLE
# WITH sex + tag-type covariates on adult survival and recapture.
#
# Three data sources sharing common demographic parameters:
#   (1) ABUNDANCE   : volunteer wetland counts, 2014-2025 (relative index)
#   (2) SURVIVAL    : adult CJS encounter histories, 2010-2026, + individual
#                     covariates (sex: F/M/Unknown; tag: none/geolocator/GPS)
#                     (the 2025-2026 interval falls outside the count series
#                     and is pooled into its own baseline, phi_post -- see 1d)
#   (3) PRODUCTIVITY: per-pair fecundity (chicks/pair), 2018-2026
#
# Survival structure:
#   logit(phi[i,t]) = lphi_base[t] + bSex_phi[sex_i] + bTag_phi[tag_i]
#   logit(p[i])     = lp_base      + bTag_p[tag_i]
#   - lphi_base[t] = time-varying UNTAGGED, FEMALE-reference survival.
#     phi_ad[t] = ilogit(lphi_base[t]) IS the population rate driving the IPM
#     process model (population is overwhelmingly untagged).
#   - Female (sex=1) and untagged (tag=1) are the reference levels (effect = 0).
#   - Tag effects constant in time (data too sparse for tag x year).
#
# Relative-scale abundance: interpret lambda / trajectory shape, not absolute N.
# Life history: breeds at age 1, two age classes (juv, adult), immigration term.
# =============================================================================

library(nimble)
library(coda)
set.seed(20260521)

## ===========================================================================
## SECTION 1 -- DATA PREP
## ===========================================================================

# ---- 1a. Counts ------------------------------------------------------------
dc <- read.csv("./Arin LEYE R code/data_2025.csv")
dc$yr_idx <- as.integer(factor(dc$year))
dc$wk_idx <- as.integer(factor(dc$week))
dc$wet_idx <- as.integer(factor(dc$wetland))
count_years <- sort(unique(dc$year)) # 2014..2025
nYearC <- length(count_years)
nWeek <- length(unique(dc$wk_idx))
nWet <- length(unique(dc$wet_idx))
nCount <- nrow(dc)
log_et <- log1p(dc$et)
log_et <- log_et - mean(log_et)
log_ed <- log1p(dc$ed)
log_ed <- log_ed - mean(log_ed)

# ---- 1b. Survival (CJS) + covariates, modern era 2010-2026 -----------------
prep_cjs_cov <- function(
  eh_path = "data/LEYE_95_26_EH.csv",
  cov_path = "data/LEYE_95_26_covs.csv"
) {
  eh <- read.csv(eh_path, check.names = FALSE)
  cov <- read.csv(cov_path, check.names = FALSE)
  nm <- names(eh)
  nm[nm == "yr_3910"] <- "yr_2010"
  nm[nm == "yr_3911"] <- "yr_2011"
  names(eh) <- nm
  stopifnot(all(eh$Band_Number == cov$Band_Number)) # required lockstep order

  yr_of <- function(c) {
    s <- sub("yr_", "", c)
    n <- as.integer(s)
    if (n >= 93 & n <= 99) {
      return(1900 + n)
    }
    if (nchar(s) == 4) {
      return(n)
    }
    2000 + n
  }
  ycols <- grep("^yr_", names(eh), value = TRUE)
  ymap <- setNames(ycols, vapply(ycols, yr_of, numeric(1)))
  years <- 2010:2026
  M <- as.matrix(eh[, ymap[as.character(years)]])
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

  # Sex==3 is the existing "Unknown" category; some newly banded individuals
  # have no Sex recorded yet (NA) rather than being coded Unknown -- treat
  # those the same way so sex[] is always a valid bSex_phi/index (1/2/3).
  sex <- as.integer(cov$Sex)
  sex[is.na(sex)] <- 3L

  list(
    y = M,
    f = f,
    nind = nrow(M),
    nocc = ncol(M),
    years = years,
    sex = sex,
    tag = as.integer(cov$Tag_type)
  )
}
cjs <- prep_cjs_cov()
cjs_years <- cjs$years
cjs_int_start <- cjs_years[-cjs$nocc] # interval start years, length nocc-1
nSex <- 3L
nTag <- 3L

# ---- 1c. Productivity (fecundity), 2018-2026 -------------------------------
dp <- read.csv("data/AnnualSucc_ProdData_95_26.csv")
names(dp) <- trimws(names(dp))
dp <- dp[!is.na(dp$SpecChickNumber), c("Year", "SpecChickNumber")]
prod_years <- sort(unique(dp$Year))
dp$py_idx <- match(dp$Year, prod_years)
nProd <- nrow(dp)
nYearP <- length(prod_years)

# ---- 1d. Master year axis & cross-dataset mappings -------------------------
ipm_years <- count_years # 2014..2025
T <- length(ipm_years)
ipm_trans_start <- ipm_years[-T] # 2014..2024 (length T-1)

# Map each CJS interval to the combined baseline-survival vector:
# lphi_vec[1..T-1] = lphi_base (untagged/female baseline, by IPM transition),
# lphi_vec[T]      = lphi_pre  (pooled pre-2014 baseline),
# lphi_vec[T+1]    = lphi_post (pooled baseline for CJS intervals starting
#                     after the count series, i.e. the 2025-2026 interval --
#                     counts haven't been updated past 2025 yet).
overlap_start <- 2014
overlap_end <- max(ipm_trans_start) # 2024
PHI_PRE_SLOT <- T
PHI_POST_SLOT <- T + 1
cjs_phi_lookup <- ifelse(
  cjs_int_start < overlap_start,
  PHI_PRE_SLOT,
  ifelse(
    cjs_int_start > overlap_end,
    PHI_POST_SLOT,
    match(cjs_int_start, ipm_trans_start)
  )
)
stopifnot(length(cjs_phi_lookup) == cjs$nocc - 1)
stopifnot(all(cjs_phi_lookup >= 1 & cjs_phi_lookup <= PHI_POST_SLOT))

# Fecundity years -> IPM year index. Years beyond the count series (2026) map
# to slot T+1, a pooled/orphan f_rate not used by the recruitment process
# model (which only ever indexes f_rate[1..T-1]) -- mirrors phi_post above.
F_POST_SLOT <- T + 1
prod_ipm_idx <- ifelse(
  prod_years <= max(ipm_years),
  match(prod_years, ipm_years), # 2018..2025 -> 5..12
  F_POST_SLOT # 2026 -> 13
)

# latent-alive inits for CJS
z_init <- function(y, f) {
  zi <- matrix(NA, nrow(y), ncol(y))
  for (i in 1:nrow(y)) {
    if (f[i] < ncol(y)) zi[i, (f[i] + 1):ncol(y)] <- 1
  }
  zi
}

## ===========================================================================
## SECTION 2 -- THE IPM
## ===========================================================================
ipmCode <- nimbleCode({
  ## ----- 2a. VITAL-RATE PRIORS ---------------------------------------------

  # --- Adult survival baseline: untagged, female-reference, time-varying ---
  mu_phi ~ dnorm(0, sd = 1.5) # mean baseline survival, logit scale
  sd_phi ~ dexp(2)
  for (t in 1:(T - 1)) {
    lphi_base[t] <- mu_phi + eps_phi[t]
    eps_phi[t] ~ dnorm(0, sd = sd_phi)
    phi_ad[t] <- ilogit(lphi_base[t]) # <-- IPM-shared POPULATION rate
    phi_juv[t] <- rho * phi_ad[t] # juvenile = fraction of adult
  }
  phi_pre ~ dunif(0, 1) # pooled pre-2014 baseline survival
  lphi_pre <- logit(phi_pre)
  phi_post ~ dunif(0, 1) # pooled baseline survival, 2025-2026 (beyond count series)
  lphi_post <- logit(phi_post)

  # --- Survival covariate effects (corner constraints: F & untagged = ref) ---
  bSex_phi[1] <- 0
  bTag_phi[1] <- 0
  for (s in 2:nSex) {
    bSex_phi[s] ~ dnorm(0, sd = 1.5)
  }
  for (g in 2:nTag) {
    bTag_phi[g] ~ dnorm(0, sd = 1.5)
  }

  # --- Recapture: constant baseline + tag effect ---
  lp_base ~ dnorm(0, sd = 1.5)
  bTag_p[1] <- 0
  for (g in 2:nTag) {
    bTag_p[g] ~ dnorm(0, sd = 1.5)
  }

  # --- Juvenile survival ratio ---
  rho ~ dbeta(2, 2)

  # --- Fecundity (chicks/pair), annual, partially pooled ---
  # f_rate[1..T] drive the recruitment process model (2c); f_rate[T+1] is the
  # pooled 2026 rate, orphaned from the process model until counts catch up.
  mu_f ~ dnorm(0, sd = 2)
  sd_f ~ dexp(1)
  for (t in 1:(T + 1)) {
    log(f_rate[t]) <- mu_f + eps_f[t]
    eps_f[t] ~ dnorm(0, sd = sd_f)
  }

  # --- Immigration ---
  omega ~ dexp(1)

  ## ----- 2b. PROCESS MODEL (population-wide, age-structured) ---------------
  N_ad[1] ~ T(dnorm(0, sd = 50), 0, )
  N_juv[1] ~ T(dnorm(0, sd = 50), 0, )
  N_tot[1] <- N_ad[1] + N_juv[1]

  for (t in 2:T) {
    R_mean[t] <- N_ad[t - 1] * (f_rate[t - 1] / 2) * phi_juv[t - 1] # recruits
    N_juv[t] ~ dpois(R_mean[t])
    S_mean[t] <- (N_ad[t - 1] + N_juv[t - 1]) * phi_ad[t - 1] # local survivors
    N_surv[t] ~ dpois(S_mean[t])
    Imm[t] ~ dpois(N_tot[t - 1] * omega) # immigrants -> adults
    N_ad[t] <- N_surv[t] + Imm[t]
    N_tot[t] <- N_ad[t] + N_juv[t]
  }
  for (t in 1:(T - 1)) {
    lambda[t] <- N_tot[t + 1] / N_tot[t]
  }

  ## ----- 2c. OBSERVATION MODEL 1: COUNTS (relative abundance) --------------
  logB0 ~ dnorm(0, sd = 5)
  b_et ~ dnorm(0, sd = 5)
  b_ed ~ dnorm(0, sd = 5)
  weekFE[1] <- 0
  for (w in 2:nWeek) {
    weekFE[w] ~ dnorm(0, sd = 5)
  }
  sd_site ~ dexp(1)
  for (j in 1:nWet) {
    z_site[j] ~ dnorm(0, sd = 1)
    site_re[j] <- z_site[j] * sd_site
  }
  sd_od ~ dexp(1)
  for (i in 1:nCount) {
    od[i] ~ dnorm(0, sd = sd_od)
    log(lam[i]) <- logB0 +
      log(N_tot[cyr[i]]) +
      weekFE[cwk[i]] +
      b_et * log_et[i] +
      b_ed * log_ed[i] +
      site_re[cwet[i]] +
      od[i]
    count[i] ~ dpois(lam[i])
  }

  ## ----- 2d. OBSERVATION MODEL 2: CJS SURVIVAL with covariates -------------
  # combined baseline (logit): slots 1..T-1 = lphi_base, slot T = lphi_pre,
  # slot T+1 = lphi_post
  for (t in 1:(T - 1)) {
    lphi_vec[t] <- lphi_base[t]
  }
  lphi_vec[T] <- lphi_pre
  lphi_vec[T + 1] <- lphi_post

  for (i in 1:nind) {
    logit(p_ind[i]) <- lp_base + bTag_p[tag[i]] # individual recapture prob
    zc[i, fcap[i]] <- 1
    for (t in (fcap[i] + 1):nocc) {
      logit(phi_ind[i, t - 1]) <- lphi_vec[cjsLookup[t - 1]] +
        bSex_phi[sex[i]] +
        bTag_phi[tag[i]]
      zc[i, t] ~ dbern(phi_ind[i, t - 1] * zc[i, t - 1])
      ych[i, t] ~ dbern(p_ind[i] * zc[i, t])
    }
  }

  ## ----- 2e. OBSERVATION MODEL 3: FECUNDITY --------------------------------
  for (m in 1:nProd) {
    chicks[m] ~ dpois(f_rate[pyr[m]])
  }
})

## ===========================================================================
## SECTION 3 -- CONSTANTS, DATA, INITS, MONITORS
## ===========================================================================
constants <- list(
  T = T,
  nWeek = nWeek,
  nWet = nWet,
  nCount = nCount,
  cyr = dc$yr_idx,
  cwk = dc$wk_idx,
  cwet = dc$wet_idx,
  nocc = cjs$nocc,
  nind = cjs$nind,
  fcap = cjs$f,
  cjsLookup = cjs_phi_lookup,
  nSex = nSex,
  nTag = nTag,
  sex = cjs$sex,
  tag = cjs$tag,
  nProd = nProd,
  pyr = prod_ipm_idx[dp$py_idx]
)
data <- list(
  count = dc$count,
  ych = cjs$y,
  chicks = dp$SpecChickNumber,
  log_et = log_et,
  log_ed = log_ed
)
inits <- function() {
  list(
    mu_phi = qlogis(0.6),
    sd_phi = 0.3,
    eps_phi = rnorm(T - 1, 0, 0.2),
    phi_pre = 0.5,
    phi_post = 0.5,
    rho = 0.4,
    bSex_phi = c(NA, 0, 0),
    bTag_phi = c(NA, 0, 0),
    lp_base = qlogis(0.5),
    bTag_p = c(NA, 0, 0),
    mu_f = log(2),
    sd_f = 0.3,
    eps_f = rnorm(T + 1, 0, 0.2),
    omega = 0.1,
    N_ad = c(40, rep(NA, T - 1)),
    N_surv = c(NA, rep(38, T - 1)),
    N_juv = rep(20, T),
    Imm = c(NA, rep(2, T - 1)),
    logB0 = -3,
    b_et = 0,
    b_ed = 0,
    weekFE = c(NA, rnorm(nWeek - 1, 0, 0.3)),
    sd_site = 0.5,
    z_site = rnorm(nWet, 0, 1),
    sd_od = 0.5,
    od = rnorm(nCount, 0, 0.3),
    zc = z_init(cjs$y, cjs$f)
  )
}
monitors <- c(
  "phi_ad",
  "phi_juv",
  "phi_pre",
  "phi_post",
  "rho",
  "bSex_phi",
  "bTag_phi",
  "bTag_p",
  "lp_base",
  "f_rate",
  "mu_f",
  "omega",
  "Imm",
  "R_mean",
  "S_mean",
  "N_surv",
  "N_tot",
  "N_ad",
  "N_juv",
  "lambda",
  "logB0",
  "sd_od",
  "sd_site"
)

## ===========================================================================
## SECTION 4 -- BUILD / RUN
## ===========================================================================
# When sourced by the parallel wrapper, SKIP_RUN is set TRUE so that only the
# model objects (ipmCode, constants, data, inits, monitors) are defined and the
# serial fit/summary below is skipped.
if (!exists("SKIP_RUN")) {
  SKIP_RUN <- FALSE
}
if (!SKIP_RUN) {
  Rmodel <- nimbleModel(
    ipmCode,
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
    niter = 10,
    nburnin = 10,
    thin = 1,
    nchains = 3,
    inits = inits,
    samplesAsCodaMCMC = TRUE,
    setSeed = 1:3
  )

  ## ===========================================================================
  ## SECTION 5 -- SUMMARY
  ## ===========================================================================
  library(MCMCvis)
  MCMCsummary(samples)

  mat <- as.matrix(samples)
  q <- function(x) {
    c(
      mean = mean(x),
      lwr = quantile(x, .025, names = FALSE),
      upr = quantile(x, .975, names = FALSE)
    )
  }

  cat("\n--- Convergence (Rhat, key params; want < 1.05) ---\n")
  keyp <- intersect(
    c(
      "phi_pre",
      "phi_post",
      "rho",
      "mu_f",
      "omega",
      "logB0",
      "lp_base",
      "bTag_phi[2]",
      "bTag_phi[3]",
      "bSex_phi[2]",
      "bSex_phi[3]",
      "bTag_p[2]",
      "bTag_p[3]"
    ),
    colnames(mat)
  )
  print(round(gelman.diag(samples[, keyp], multivariate = FALSE)$psrf, 3))

  cat(
    "\n--- Baseline adult survival phi_ad (untagged, female) by transition ---\n"
  )
  pa <- t(apply(mat[, paste0("phi_ad[", 1:(T - 1), "]")], 2, q))
  print(
    data.frame(
      transition = paste0(ipm_years[-T], "->", ipm_years[-1]),
      round(pa, 3)
    ),
    row.names = FALSE
  )
  cat(sprintf(
    "pooled pre-2014 (phi_pre): %.3f | pooled 2025->2026 (phi_post): %.3f\n",
    mean(mat[, "phi_pre"]),
    mean(mat[, "phi_post"])
  ))

  cat("\n--- Survival covariate effects (logit scale) and odds ratios ---\n")
  eff <- rbind(
    Sex_M_vs_F = q(mat[, "bSex_phi[2]"]),
    Sex_Unk_vs_F = q(mat[, "bSex_phi[3]"]),
    Tag_geo_vs_none = q(mat[, "bTag_phi[2]"]),
    Tag_GPS_vs_none = q(mat[, "bTag_phi[3]"])
  )
  print(round(eff, 3))
  cat("Odds ratios (exp):\n")
  cat(sprintf(
    "  geolocator vs untagged: %.2f (%.2f, %.2f)\n",
    exp(mean(mat[, "bTag_phi[2]"])),
    exp(quantile(mat[, "bTag_phi[2]"], .025)),
    exp(quantile(mat[, "bTag_phi[2]"], .975))
  ))
  cat(sprintf(
    "  GPS vs untagged:        %.2f (%.2f, %.2f)\n",
    exp(mean(mat[, "bTag_phi[3]"])),
    exp(quantile(mat[, "bTag_phi[3]"], .025)),
    exp(quantile(mat[, "bTag_phi[3]"], .975))
  ))
  cat(sprintf(
    "P(geolocator reduces survival): %.3f | P(GPS reduces survival): %.3f\n",
    mean(mat[, "bTag_phi[2]"] < 0),
    mean(mat[, "bTag_phi[3]"] < 0)
  ))

  cat("\n--- Recapture: baseline + tag effects on p ---\n")
  cat(sprintf("baseline p (untagged): %.3f\n", plogis(mean(mat[, "lp_base"]))))
  cat(sprintf(
    "  tag effect on p (geo): OR %.2f | (GPS): OR %.2f\n",
    exp(mean(mat[, "bTag_p[2]"])),
    exp(mean(mat[, "bTag_p[3]"]))
  ))

  cat("\n--- Fecundity f_rate (chicks/pair) by year ---\n")
  fr <- t(apply(mat[, paste0("f_rate[", 1:T, "]")], 2, q))
  print(data.frame(year = ipm_years, round(fr, 3)), row.names = FALSE)
  cat(sprintf(
    "pooled 2026 (beyond count series): %.3f\n",
    mean(mat[, paste0("f_rate[", T + 1, "]")])
  ))

  cat("\n--- Population growth lambda ---\n")
  lam <- t(apply(mat[, paste0("lambda[", 1:(T - 1), "]")], 2, q))
  print(
    data.frame(
      transition = paste0(ipm_years[-T], "->", ipm_years[-1]),
      round(lam, 3)
    ),
    row.names = FALSE
  )
  geolam <- exp(mean(log(rowMeans(mat[, paste0("lambda[", 1:(T - 1), "]")]))))
  cat(sprintf(
    "Approx geometric-mean lambda: %.3f (%.1f%%/yr)\n",
    geolam,
    100 * (geolam - 1)
  ))

  cat("\n--- Immigration omega ---\n")
  print(round(q(mat[, "omega"]), 3))
  cat(sprintf("P(omega > 0.01): %.3f\n", mean(mat[, "omega"] > 0.01)))

  # ---- group survival helper (sex x tag) -------------------------------------
  group_survival <- function(mat, year_idx = 1) {
    base <- qlogis(mat[, paste0("phi_ad[", year_idx, "]")])
    sexlab <- c("F", "M", "Unk")
    taglab <- c("none", "geo", "GPS")
    bS <- cbind(0, mat[, "bSex_phi[2]"], mat[, "bSex_phi[3]"])
    bT <- cbind(0, mat[, "bTag_phi[2]"], mat[, "bTag_phi[3]"])
    out <- list()
    for (s in 1:3) {
      for (g in 1:3) {
        v <- plogis(base + bS[, s] + bT[, g])
        out[[paste(sexlab[s], taglab[g], sep = "_")]] <-
          c(
            mean = mean(v),
            lwr = quantile(v, .025, names = FALSE),
            upr = quantile(v, .975, names = FALSE)
          )
      }
    }
    round(do.call(rbind, out), 3)
  }
  cat(
    "\n--- Survival by sex x tag (transition 1; * GPS x Unknown has no data) ---\n"
  )
  print(group_survival(mat, 1))

  # ---- abundance trajectory plot ---------------------------------------------
  Nt <- t(apply(mat[, paste0("N_tot[", 1:T, "]")], 2, q))
  #png("leye_ipm_abundance.png", width=1600, height=1000, res=200)
  par(mar = c(4.5, 4.5, 2.5, 1))
  plot(
    ipm_years,
    Nt[, 1],
    type = "n",
    ylim = range(c(Nt[, 2], Nt[, 3])),
    xlab = "Year",
    ylab = "Latent total abundance (relative scale)",
    main = "LEYE IPM: estimated population trajectory"
  )
  polygon(
    c(ipm_years, rev(ipm_years)),
    c(Nt[, 2], rev(Nt[, 3])),
    col = adjustcolor("firebrick", .18),
    border = NA
  )
  lines(ipm_years, Nt[, 1], col = "firebrick", lwd = 2)
  points(ipm_years, Nt[, 1], pch = 16, col = "firebrick")
  #dev.off()
  cat("\nWrote: leye_ipm_abundance.png\n")

  # =============================================================================
  # Forest plot: sex + tag covariate effects from the LEYE IPM
  # Top panel : effects on adult SURVIVAL (phi)   -- odds-ratio scale
  # Bottom    : effect of tag on RECAPTURE (p)     -- odds-ratio scale
  # Reference levels (female, untagged) shown as a dashed line at OR = 1.
  # Assumes `samples` (coda) is in the workspace.
  # =============================================================================

  mat <- as.matrix(samples)
  orq <- function(col) {
    # odds ratio summary for a logit coef
    x <- exp(mat[, col])
    c(
      med = median(x),
      lwr = quantile(x, .025, names = FALSE),
      upr = quantile(x, .975, names = FALSE),
      pneg = mean(mat[, col] < 0)
    ) # P(effect reduces the rate)
  }

  # ---- survival effects (vs reference: female, untagged) ----
  phi_rows <- rbind(
    "Male (vs female)" = orq("bSex_phi[2]"),
    "Unknown sex (vs female)" = orq("bSex_phi[3]"),
    "Geolocator (vs untagged)" = orq("bTag_phi[2]"),
    "GPS (vs untagged)" = orq("bTag_phi[3]")
  )
  # ---- recapture effects (tag, vs untagged) ----
  p_rows <- rbind(
    "Geolocator (vs untagged)" = orq("bTag_p[2]"),
    "GPS (vs untagged)" = orq("bTag_p[3]")
  )

  #png("leye_covariate_forest.png", width = 1700, height = 1150, res = 200)
  layout(matrix(1:2, nrow = 2), heights = c(4, 2.6))
  par(mar = c(3, 12, 3, 6), xaxs = "i")

  draw_forest <- function(tab, title, xlim, col_pts) {
    n <- nrow(tab)
    ypos <- n:1
    plot(
      NA,
      xlim = xlim,
      ylim = c(0.5, n + 0.5),
      log = "x",
      xlab = "",
      ylab = "",
      yaxt = "n",
      main = title
    )
    abline(v = 1, lty = 2, col = "grey50") # OR = 1 (no effect)
    segments(tab[, "lwr"], ypos, tab[, "upr"], ypos, col = col_pts, lwd = 2)
    points(tab[, "med"], ypos, pch = 16, col = col_pts, cex = 1.4)
    axis(2, at = ypos, labels = rownames(tab), las = 1, cex.axis = 0.9)
    # annotate P(reduces rate) on the right margin
    for (k in 1:n) {
      mtext(
        sprintf("P(<1)=%.2f", tab[k, "pneg"]),
        side = 4,
        at = ypos[k],
        las = 1,
        line = 0.3,
        cex = 0.7,
        col = "grey30"
      )
    }
  }

  # shared x-limits so panels are comparable; pad around observed range
  allor <- rbind(phi_rows, p_rows)
  xl <- range(c(allor[, "lwr"], allor[, "upr"]))
  xl <- c(min(0.3, xl[1] * 0.9), max(3, xl[2] * 1.1))

  draw_forest(phi_rows, "Effect on adult survival (odds ratio)", xl, "#1f4e79")
  mtext(
    "OR < 1 = lower survival than reference",
    side = 1,
    line = 2,
    cex = 0.7,
    col = "grey30"
  )
  draw_forest(
    p_rows,
    "Effect on recapture probability (odds ratio)",
    xl,
    "#7030a0"
  )
  mtext(
    "Odds ratio (log scale); points = posterior median, bars = 95% CrI",
    side = 1,
    line = 2,
    cex = 0.7,
    col = "grey30"
  )
  #dev.off()
  cat("Wrote: leye_covariate_forest.png\n\n")

  cat("Survival effects (odds-ratio scale):\n")
  print(round(phi_rows, 3))
  cat("\nRecapture effects (odds-ratio scale):\n")
  print(round(p_rows, 3))

  # =============================================================================
  # Figure: adult vs juvenile annual survival from the LEYE IPM (covariate model)
  #
  # phi_ad[t]  = baseline adult survival = untagged, FEMALE-reference birds.
  #              (This is the population rate driving the process model.)
  # phi_juv[t] = rho * phi_ad[t]  (juvenile survival; rho weakly identified.)
  #
  # Tagged/male survival differs by the fitted covariate effects; this figure
  # shows the BASELINE population rates, not tagged-bird rates. Use the forest
  # plot (leye_covariate_forest.R) for the covariate effects themselves.
  #
  # Assumes `samples` (coda) and `ipm_years` are in the workspace.
  # =============================================================================

  mat <- as.matrix(samples)
  Tm1 <- length(ipm_years) - 1
  trans_yr <- ipm_years[-length(ipm_years)] # start year of each transition

  summ <- function(stem) {
    cols <- paste0(stem, "[", 1:Tm1, "]")
    t(apply(mat[, cols, drop = FALSE], 2, function(x) {
      c(
        med = median(x),
        lwr = quantile(x, .025, names = FALSE),
        upr = quantile(x, .975, names = FALSE)
      )
    }))
  }
  ad <- summ("phi_ad")
  juv <- summ("phi_juv")

  # transitions with direct CJS data (start year present in CJS intervals)
  has_cjs <- trans_yr %in% cjs_int_start

  #png("leye_survival_adult_juvenile.png", width = 1700, height = 1200, res = 200)\
  layout(matrix(1, nrow = 1), heights = c(6))
  par(mar = c(5.5, 4.5, 3, 1))

  dodge <- 0.05
  xa <- trans_yr - dodge
  xj <- trans_yr + dodge
  col_ad <- "#1f4e79"
  col_juv <- "#c0504d"

  plot(
    NA,
    xlim = range(trans_yr) + c(-0.5, 0.5),
    ylim = c(0, 1),
    xlab = "",
    ylab = "Annual apparent survival",
    main = "LEYE adult vs juvenile survival (baseline: untagged, female)",
    xaxt = "n"
  )
  axis(1, at = trans_yr, labels = paste0(trans_yr + 1), cex.axis = 1, las = 2)
  abline(h = seq(0, 1, .2), col = "grey92")

  # adult
  arrows(
    xa,
    ad[, "lwr"],
    xa,
    ad[, "upr"],
    length = .025,
    angle = 90,
    code = 3,
    col = col_ad,
    lwd = 1.5
  )
  points(xa, ad[, "med"], pch = 16, col = col_ad, cex = 1.2)
  lines(xa, ad[, "med"], col = col_ad, lty = 2)
  # juvenile
  arrows(
    xj,
    juv[, "lwr"],
    xj,
    juv[, "upr"],
    length = .025,
    angle = 90,
    code = 3,
    col = col_juv,
    lwd = 1.5
  )
  points(xj, juv[, "med"], pch = 17, col = col_juv, cex = 1.2)
  lines(xj, juv[, "med"], col = col_juv, lty = 2)

  # mark transitions lacking direct CJS data (informed only via process model)
  no_dat <- trans_yr[!has_cjs]
  if (length(no_dat)) {
    points(no_dat, rep(0, length(no_dat)), pch = 8, col = "grey40", cex = 0.8)
  }

  legend(
    "topright",
    bty = "n",
    cex = 0.9,
    legend = c("Adult (phi_ad)", "Juvenile (phi_juv = rho x phi_ad)"),
    col = c(col_ad, col_juv),
    pch = c(16, 17),
    lty = 2,
    lwd = 1.5
  )
  mtext(
    "Points = posterior median; bars = 95% CrI.   * (grey) = transition without direct CJS data.",
    side = 1,
    line = 3.6,
    cex = 0.62,
    adj = 0
  )
  mtext(
    "Juvenile band inherits adult variation (shared rho); do not read year-to-year juvenile trend as independent.",
    side = 1,
    line = 4.2,
    cex = 0.62,
    adj = 0
  )
  #dev.off()
  cat("Wrote: leye_survival_adult_juvenile.png\n")

  cat("\nAdult survival (baseline):\n")
  print(
    data.frame(transition = paste0(trans_yr, "-", trans_yr + 1), round(ad, 3)),
    row.names = FALSE
  )
  cat("\nJuvenile survival:\n")
  print(
    data.frame(transition = paste0(trans_yr, "-", trans_yr + 1), round(juv, 3)),
    row.names = FALSE
  )

  # =============================================================================
  # Figure: annual adult vs juvenile survival from the LEYE IPM
  # Plots posterior median + 95% CrI for phi_ad[t] and phi_juv[t] by transition year.
  # Assumes `samples` (coda mcmc.list) and `ipm_years` are in the workspace.
  # =============================================================================

  mat <- as.matrix(samples)
  Tm1 <- length(ipm_years) - 1 # number of survival transitions
  trans_yr <- ipm_years[-length(ipm_years)] # start year of each transition

  # pull posterior summaries for a parameter vector par[1..Tm1]
  summ <- function(stem) {
    cols <- paste0(stem, "[", 1:Tm1, "]")
    t(apply(mat[, cols, drop = FALSE], 2, function(x) {
      c(
        med = median(x),
        lwr = quantile(x, .025, names = FALSE),
        upr = quantile(x, .975, names = FALSE)
      )
    }))
  }
  ad <- summ("phi_ad")
  juv <- summ("phi_juv")

  # CJS overlap years carry direct survival data (2014-2022 starts); flag for the reader
  has_cjs <- trans_yr %in% (cjs_int_start) # cjs_int_start from data prep

  #png("leye_survival_by_year.png", width = 1700, height = 1050, res = 200)
  par(mar = c(4.5, 4.5, 2.5, 1))

  dodge <- 0.05 # offset so CIs don't overlap
  xa <- trans_yr - dodge
  xj <- trans_yr + dodge
  col_ad <- "#1f4e79"
  col_juv <- "#c0504d"

  plot(
    NA,
    xlim = range(trans_yr) + c(-0.5, 0.5),
    ylim = c(0, 1),
    xlab = "Transition (survival from year t to t+1)",
    ylab = "Annual apparent survival",
    main = "LEYE adult vs juvenile survival (IPM posterior)",
    xaxt = "n"
  )
  axis(1, at = trans_yr, labels = paste0(trans_yr + 1), cex.axis = 1, las = 2)
  abline(h = seq(0, 1, .2), col = "grey92")

  # adult: CI bars + points
  arrows(
    xa,
    ad[, "lwr"],
    xa,
    ad[, "upr"],
    length = .025,
    angle = 90,
    code = 3,
    col = col_ad,
    lwd = 1.5
  )
  points(xa, ad[, "med"], pch = 16, col = col_ad, cex = 1.2)
  lines(xa, ad[, "med"], col = col_ad, lwd = 1, lty = 2)

  # juvenile: CI bars + points
  arrows(
    xj,
    juv[, "lwr"],
    xj,
    juv[, "upr"],
    length = .025,
    angle = 90,
    code = 3,
    col = col_juv,
    lwd = 1.5
  )
  points(xj, juv[, "med"], pch = 17, col = col_juv, cex = 1.2)
  lines(xj, juv[, "med"], col = col_juv, lwd = 1, lty = 2)

  # mark transitions lacking direct CJS data (juvenile always indirect; adult here too)
  no_dat <- trans_yr[!has_cjs]
  if (length(no_dat)) {
    mtext("*", side = 1, at = no_dat, line = -0.5, col = "grey40", cex = 1.5)
  }

  legend(
    "topright",
    bty = "n",
    cex = .9,
    legend = c(
      "Adult survival (phi_ad)",
      "Juvenile survival (phi_juv = rho x phi_ad)"
    ),
    col = c(col_ad, col_juv),
    pch = c(16, 17),
    lty = 2,
    lwd = 1.5
  )
  mtext(
    "Points = posterior median; bars = 95% CrI\n* = transition without direct CJS data (informed via process model).",
    side = 3,
    line = -1.6,
    cex = 0.62,
    adj = 0.02
  )
  #dev.off()
  cat("Wrote: leye_survival_by_year.png\n")

  # also print the table behind the figure
  cat("\nAdult survival:\n")
  print(
    data.frame(transition = paste0(trans_yr, "-", trans_yr + 1), round(ad, 3)),
    row.names = FALSE
  )
  cat("\nJuvenile survival:\n")
  print(
    data.frame(transition = paste0(trans_yr, "-", trans_yr + 1), round(juv, 3)),
    row.names = FALSE
  )

  # =============================================================================
  # Retrospective Analyses --------------------------------------------------
  # =============================================================================
  nimble_to_simslist <- function(samples) {
    # Accepts: an mcmc.list (runMCMC(..., samplesAsCodaMCMC = TRUE)),
    #          a plain list of per-chain matrices (samplesAsCodaMCMC = FALSE, multi-chain),
    #          or a single matrix/mcmc object (one chain).
    if (inherits(samples, "mcmc.list")) {
      mat <- as.matrix(samples) # rbinds chains, keeps varnames as colnames
    } else if (is.list(samples) && !is.data.frame(samples)) {
      mat <- do.call(rbind, lapply(samples, as.matrix))
    } else if (is.matrix(samples) || inherits(samples, "mcmc")) {
      mat <- as.matrix(samples)
    } else {
      stop("`samples` must be an mcmc.list, a list of matrices, or a matrix.")
    }

    nms <- colnames(mat)
    if (is.null(nms)) {
      stop("Columns must be named, e.g. 'N[1,2]'.")
    }
    n <- nrow(mat) # total post-burnin/thin draws across chains

    has_idx <- grepl("\\[", nms)
    base <- ifelse(has_idx, sub("\\[.*$", "", nms), nms)

    # Parse each column's integer index vector (empty for scalars)
    idx <- lapply(seq_along(nms), function(i) {
      if (!has_idx[i]) {
        return(integer(0))
      }
      inside <- sub(".*\\[(.*)\\]$", "\\1", nms[i])
      as.integer(strsplit(inside, ",", fixed = TRUE)[[1]])
    })

    out <- list()
    for (b in unique(base)) {
      cols <- which(base == b)

      if (length(idx[[cols[1]]]) == 0L) {
        # scalar parameter -> length-n vector
        out[[b]] <- mat[, cols[1]]
        next
      }

      idx_mat <- do.call(rbind, idx[cols]) # (n_cols x n_dims) parsed indices
      n_dims <- ncol(idx_mat)
      dims <- apply(idx_mat, 2, max) # extent per parameter dimension
      arr <- array(NA_real_, dim = c(n, dims)) # iteration dim first

      # Fill each column into its exact array slice via matrix indexing
      for (k in seq_along(cols)) {
        coord <- cbind(
          seq_len(n),
          matrix(idx_mat[k, ], nrow = n, ncol = n_dims, byrow = TRUE)
        )
        arr[coord] <- mat[, cols[k]]
      }
      out[[b]] <- arr
    }
    out
  } # custom function to convert nimble output to sims.list object

  draws <- nimble_to_simslist(samples)
  n.draws <- nrow(mat) # Determine number of MCMC draws
  corr <- matrix(NA, ncol = 3, nrow = n.draws) # Create object to hold results
  for (s in 1:n.draws) {
    # Loop over all MCMC draws and get correlations
    corr[s, 1] <- cor(draws$phi_juv[s, ], draws$lambda[s, ])
    corr[s, 2] <- cor(draws$phi_ad[s, ], draws$lambda[s, ])
    corr[s, 3] <- cor(draws$f_rate[s, 1:11], draws$lambda[s, ])
  }

  # Calculate posterior summaries for the correlation coefficients
  # Posterior means
  apply(corr, 2, mean)
  # 95% credible intervals
  cri <- function(x) quantile(x, c(0.025, 0.975))
  apply(corr, 2, cri)

  # ~~~~ code for Figure 9.2 ~~~~
  cri <- function(x) quantile(x, c(0.025, 0.975))
  cri.lambda <- apply(draws$lambda, 2, cri)

  op <- par(mfrow = c(2, 3))
  layout(
    matrix(1:3, 1, 3, byrow = TRUE),
    widths = c(3.05, 3, 3),
    heights = rep(3.5, 3),
    TRUE
  )

  cri.rate <- apply(draws$phi_juv, 2, cri)
  plot(
    NA,
    ylim = range(cri.lambda),
    xlim = range(cri.rate),
    ylab = expression('Population growth rate (' * lambda * ')'),
    xlab = expression('Juvenile survival (' * phi[italic(j)] * ')'),
    axes = FALSE
  )
  axis(1)
  axis(1, at = c(0.075, 0.125, 0.175, 0.225), labels = NA, tcl = -0.25)
  axis(2, las = 1)
  segments(
    colMeans(draws$phi_juv),
    cri.lambda[1, ],
    colMeans(draws$phi_juv),
    cri.lambda[2, ],
    col = "grey60"
  )
  segments(
    cri.rate[1, ],
    colMeans(draws$lambda),
    cri.rate[2, ],
    colMeans(draws$lambda),
    col = "grey60"
  )
  points(y = colMeans(draws$lambda), x = colMeans(draws$phi_juv), pch = 16)

  cri.rate <- apply(draws$phi_ad, 2, cri)
  plot(
    NA,
    ylim = range(cri.lambda),
    xlim = range(cri.rate),
    ylab = NA,
    xlab = expression('Adult survival (' * phi[italic(a)] * ')'),
    axes = FALSE
  )
  axis(1)
  axis(2, las = 1)
  segments(
    colMeans(draws$phi_ad),
    cri.lambda[1, ],
    colMeans(draws$phi_ad),
    cri.lambda[2, ],
    col = "grey60"
  )
  segments(
    cri.rate[1, ],
    colMeans(draws$lambda),
    cri.rate[2, ],
    colMeans(draws$lambda),
    col = "grey60"
  )
  points(y = colMeans(draws$lambda), x = colMeans(draws$phi_ad), pch = 16)

  cri.rate <- apply(draws$f_rate, 2, cri)
  plot(
    NA,
    ylim = range(cri.lambda),
    xlim = range(cri.rate),
    ylab = NA,
    xlab = expression('Productivity ad (' * italic(f)[1] * ')'),
    axes = FALSE
  )
  axis(1)
  axis(2, las = 1)
  segments(
    colMeans(draws$f_rate)[1:11],
    cri.lambda[1, ],
    colMeans(draws$f_rate)[1:11],
    cri.lambda[2, ],
    col = "grey60"
  )
  segments(
    cri.rate[1, 1:11],
    colMeans(draws$lambda),
    cri.rate[2, 1:11],
    colMeans(draws$lambda),
    col = "grey60"
  )
  points(y = colMeans(draws$lambda), x = colMeans(draws$f_rate)[1:11], pch = 16)

  # LTRE --------------------------------------------------------------------
  # Calculation of growth rate sensitivities
  delta <- 0.001 # (Small) size of perturbation
  T <- 5 # Number of years with projections
  n.draws <- nrow(mat) # Number of MCMC draws
  # Define arrays to store output
  N.ref <- N.star.phi_juv <- N.star.phi_ad <- N.star.f_rate <- array(
    NA,
    dim = c(n.draws, 3, T)
  )
  N.ref[,, 1] <- N.star.phi_juv[,, 1] <- N.star.phi_ad[,, 1] <- N.star.f_rate[,,
    1
  ] <- 1
  r.annual.ref <- r.annual.star.phi_juv <- r.annual.star.phi_ad <- r.annual.star.f_rate <- matrix(
    NA,
    nrow = n.draws,
    ncol = T
  )
  lambda <- numeric()
  sens <- matrix(NA, nrow = n.draws, ncol = 5)
  # draws <- out1$sims.list # store MCMC draws in a new object to simplify calculation
  # Loop over all MCMC draws and project in time
  for (s in 1:n.draws) {
    # Loop over all MCMC draws
    for (t in 1:(T - 1)) {
      # Loop over all time steps
      # Calculate asymptotic population growth rate based on stage-specific population sizes
      N.ref[s, 1, t + 1] <- mean(draws$phi_juv[s, ]) *
        (draws$mean.f1[s] *
          (N.ref[s, 1, t] + N.ref[s, 3, t]) +
          draws$mean.f2[s] * N.ref[s, 2, t]) /
        2
      N.ref[s, 2, t + 1] <- draws$mean.phia[s] * sum(N.ref[s, , t])
      N.ref[s, 3, t + 1] <- draws$mean.omega[s] * sum(N.ref[s, , t])
      r.annual.ref[s, t] <- log(sum(N.ref[s, , t + 1])) -
        log(sum(N.ref[s, , t]))
      # Sensitivity with respect to juvenile survival
      N.star.phij[s, 1, t + 1] <- (draws$mean.phij[s] + delta) /
        2 *
        (draws$mean.f1[s] *
          (N.star.phij[s, 1, t] + N.star.phij[s, 3, t]) +
          draws$mean.f2[s] * N.star.phij[s, 2, t])
      N.star.phij[s, 2, t + 1] <- draws$mean.phia[s] * sum(N.star.phij[s, , t])
      N.star.phij[s, 3, t + 1] <- draws$mean.omega[s] * sum(N.star.phij[s, , t])
      r.annual.star.phij[s, t] <- log(sum(N.star.phij[s, , t + 1])) -
        log(sum(N.star.phij[s, , t]))
      # Sensitivity with respect to adult survival
      N.star.phia[s, 1, t + 1] <- draws$mean.phij[s] /
        2 *
        (draws$mean.f1[s] *
          (N.star.phia[s, 1, t] + N.star.phia[s, 3, t]) +
          draws$mean.f2[s] * N.star.phia[s, 2, t])
      N.star.phia[s, 2, t + 1] <- (draws$mean.phia[s] + delta) *
        sum(N.star.phia[s, , t])
      N.star.phia[s, 3, t + 1] <- draws$mean.omega[s] * sum(N.star.phia[s, , t])
      r.annual.star.phia[s, t] <- log(sum(N.star.phia[s, , t + 1])) -
        log(sum(N.star.phia[s, , t]))

      # Sensitivity with respect to immigration
      N.star.om[s, 1, t + 1] <- draws$mean.phij[s] /
        2 *
        (draws$mean.f1[s] *
          (N.star.om[s, 1, t] +
            N.star.om[s, 3, t]) +
          draws$mean.f2[s] * N.star.om[s, 2, t])
      N.star.om[s, 2, t + 1] <- draws$mean.phia[s] * sum(N.star.om[s, , t])
      N.star.om[s, 3, t + 1] <- (draws$mean.omega[s] + delta) *
        sum(N.star.om[s, , t])
      r.annual.star.om[s, t] <- log(sum(N.star.om[s, , t + 1])) -
        log(sum(N.star.om[s, , t]))
      # Sensitivity with respect to productivity 1y
      N.star.f1[s, 1, t + 1] <- draws$mean.phij[s] /
        2 *
        ((draws$mean.f1[s] + delta) *
          (N.star.f1[s, 1, t] + N.star.f1[s, 3, t]) +
          draws$mean.f2[s] * N.star.f1[s, 2, t])
      N.star.f1[s, 2, t + 1] <- draws$mean.phia[s] * sum(N.star.f1[s, , t])
      N.star.f1[s, 3, t + 1] <- draws$mean.omega[s] * sum(N.star.f1[s, , t])
      r.annual.star.f1[s, t] <- log(sum(N.star.f1[s, , t + 1])) -
        log(sum(N.star.f1[s, , t]))
      # Sensitivity with respect to productivity ad
      N.star.f2[s, 1, t + 1] <- draws$mean.phij[s] /
        2 *
        (draws$mean.f1[s] *
          (N.star.f2[s, 1, t] +
            N.star.f2[s, 3, t]) +
          (draws$mean.f2[s] + delta) * N.star.f2[s, 2, t])
      N.star.f2[s, 2, t + 1] <- draws$mean.phia[s] * sum(N.star.f2[s, , t])
      N.star.f2[s, 3, t + 1] <- draws$mean.omega[s] * sum(N.star.f2[s, , t])
      r.annual.star.f2[s, t] <- log(sum(N.star.f2[s, , t + 1])) -
        log(sum(N.star.f2[s, , t]))
    } #t
    lambda[s] <- exp(r.annual.ref[s, T - 1])
    # Growth rate sensitivities
    sens[s, 1] <- (exp(r.annual.star.phij[s, T - 1]) - lambda[s]) / delta
    sens[s, 2] <- (exp(r.annual.star.phia[s, T - 1]) - lambda[s]) / delta
    sens[s, 3] <- (exp(r.annual.star.om[s, T - 1]) - lambda[s]) / delta
    sens[s, 4] <- (exp(r.annual.star.f1[s, T - 1]) - lambda[s]) / delta
    sens[s, 5] <- (exp(r.annual.star.f2[s, T - 1]) - lambda[s]) / delta
  } #s
} # end if(!SKIP_RUN)
