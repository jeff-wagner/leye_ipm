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
#
# Life history / recruitment:
#   Chicks are absent from the breeding site until they recruit, and the 1990s
#   chick data show recruitment split roughly evenly between ages 1 and 2. So
#   the process model carries a two-year lag: recruits in year t come from the
#   cohorts fledged in t-1 (prob psi1) and t-2 (prob psi2), entering the
#   breeding class directly. There is no separate juvenile stage in N_tot,
#   because volunteer counts cannot distinguish non-breeding age-1 birds from
#   breeding adults.
#
#   psi1/psi2 are COMPOSITE survive-AND-return-locally probabilities, not
#   annual survival rates. They replace the old phi_juv = rho * phi_ad, whose
#   rho was weakly identified because recruitment and immigration (omega) both
#   add birds and the relative-scale counts cannot separate them. Priors come
#   from scripts/LEYE_juv90_module.R and are informative but updatable.
# =============================================================================

library(nimble)
library(coda)
set.seed(20260521)

## ===========================================================================
# SECTION 1 -- DATA PREP -------------------------------------------------
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
# SECTION 2 -- THE IPM ---------------------------------------------------
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

  # --- Local recruitment of fledged chicks -----------------------------------
  # psi1/psi2 = P(a fledged chick survives AND returns to breed locally at age
  # 1 / age 2). These are COMPOSITE survival x philopatry probabilities, not
  # annual survival rates. Birds that survive but settle elsewhere are not
  # counted here; immigration from other populations enters via omega.
  #
  # Priors are Beta fits to the posteriors from the 1990s chick+adult module
  # (scripts/LEYE_juv90_module.R): 139 chicks banded 1995-1998 analysed jointly
  # with 178 same-era adults, with recruitment age marginalised over.
  # Priors on pi_rec/kappa rather than psi1/psi2 directly because the former
  # are essentially uncorrelated in that posterior (r = -0.01), so independent
  # priors do not distort the joint. Informative but updatable.
  pi_rec ~ dbeta(2.268, 13.610) # total P(recruit locally): 1990s mean 0.143
  kappa ~ dbeta(6.054, 6.309) # fraction recruiting at age 1: 1990s mean 0.490
  psi1 <- pi_rec * kappa
  psi2 <- pi_rec * (1 - kappa)

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
  # Constant per-capita immigration rate. This is a TESTED assumption, not a
  # default -- do not "improve" it with a year random effect without reading
  # the following.
  #
  # A time-varying version was fitted (2026-08-05):
  #   log(omega[t]) = mu_omega + z_omega[t] * sd_omega,  sd_omega ~ dexp(2)
  # It converged cleanly (Rhat <= 1.012) and was uninformative:
  #
  #   * sd_omega posterior ~= its prior. P(sd_omega > 0.25) was 0.607 under the
  #     prior and 0.572 after fitting; prior mean 0.500 -> posterior 0.429.
  #   * Annual rates came out flat: 0.069-0.088, a 1.27x range, with every
  #     interval spanning ~0.004-0.35.
  #   * In the tLTRE the driver-explained share of Var(lambda) rose 20.9% ->
  #     28.2% while the demographic + cross terms fell by the same 7.3 points.
  #     The year effect was absorbing Poisson noise and re-presenting it as
  #     environmental variation, not detecting signal.
  #   * It also shifted the level without evidence: omega fell and pi_rec rose
  #     to compensate, because the two enter lambda additively and the counts
  #     cannot separate them.
  #
  # The cause is arithmetic, not fixable by reparameterisation: ~3 immigrants
  # arrive per year, so a rate estimated from one year of data carries a 95%
  # interval spanning roughly a 14-fold range, and the count series supplies
  # only ~10 usable growth increments in total. Resolving annual immigration
  # needs data that separates immigrants from local recruits -- chick-origin
  # returns versus newly captured unbanded adults -- not a richer model.
  omega ~ dexp(1)

  ## ----- 2b. PROCESS MODEL (breeding population, delayed recruitment) ------
  # N_ad[t] = birds present in the local BREEDING population in year t.
  # Chicks are absent until they recruit (age 1 or 2), so they are not a
  # separate observable class: recruits enter N_ad directly. Volunteer counts
  # cannot distinguish non-breeding age-1 birds from breeding adults, so
  # N_tot == N_ad rather than adding an age-1 stage to the observed total.
  #
  # Two initial years are free because recruitment carries a two-year lag.
  N_ad[1] ~ T(dnorm(0, sd = 50), 0, )
  N_ad[2] ~ T(dnorm(0, sd = 50), 0, )
  N_tot[1] <- N_ad[1]
  N_tot[2] <- N_ad[2]

  for (t in 3:T) {
    # recruits from the cohorts fledged 1 and 2 years ago
    R_mean[t] <- N_ad[t - 1] *
      (f_rate[t - 1] / 2) *
      psi1 +
      N_ad[t - 2] * (f_rate[t - 2] / 2) * psi2
    N_rec[t] ~ dpois(R_mean[t])
    S_mean[t] <- N_ad[t - 1] * phi_ad[t - 1] # local survivors
    N_surv[t] ~ dpois(S_mean[t])
    Imm[t] ~ dpois(N_ad[t - 1] * omega) # immigrants -> breeding pop
    N_ad[t] <- N_surv[t] + N_rec[t] + Imm[t]
    N_tot[t] <- N_ad[t]
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
# SECTION 3 -- CONSTANTS, DATA, INITS, MONITORS --------------------------
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
    pi_rec = 0.14,
    kappa = 0.5,
    bSex_phi = c(NA, 0, 0),
    bTag_phi = c(NA, 0, 0),
    lp_base = qlogis(0.5),
    bTag_p = c(NA, 0, 0),
    mu_f = log(2),
    sd_f = 0.3,
    eps_f = rnorm(T + 1, 0, 0.2),
    omega = 0.1,
    N_ad = c(40, 40, rep(NA, T - 2)),
    N_surv = c(NA, NA, rep(30, T - 2)),
    N_rec = c(NA, NA, rep(6, T - 2)),
    Imm = c(NA, NA, rep(4, T - 2)),
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
  "phi_pre",
  "phi_post",
  "pi_rec",
  "kappa",
  "psi1",
  "psi2",
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
  "N_rec",
  "N_tot",
  "N_ad",
  "lambda",
  "logB0",
  "sd_od",
  "sd_site"
)

## ===========================================================================
# SECTION 4 -- BUILD / RUN -----------------------------------------------
## ===========================================================================
# When sourced by the parallel wrapper, SKIP_RUN is set TRUE so that only the
# model objects (ipmCode, constants, data, inits, monitors) are defined and the
# serial fit/summary below is skipped.
if (!exists("SKIP_RUN")) {
  SKIP_RUN <- FALSE
}

## ===========================================================================
## SECTION 5 -- CONSOLE SUMMARY
##
## This is the QUICK-LOOK summary, for checking a fit at the console. The
## shareable, collaborator-facing artefact is the Quarto report:
##
##     quarto render reports/leye_ipm_report.qmd
##
## which reads the cached posterior, embeds the figures, and carries the
## interpretation and caveats. This function deliberately no longer writes
## HTML -- having two HTML generators meant two things to keep in sync, and
## the one here could not embed its own figures.
##
## Defined outside the `if (!SKIP_RUN)` block so it is available when this
## script is sourced, including by the parallel wrapper, which calls it on
## the merged multi-chain `samples`.
## ===========================================================================
summarize_ipm <- function(
  samples,
  ipm_years,
  T,
  cjs_int_start,
  plot_dir = "output/plots"
) {
  library(MCMCvis)
  library(gt)

  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }

  # Render a summary vector/matrix/data.frame as a table. Matrices with row
  # names (e.g. from rbind()/apply()) get those names as a leading column;
  # plain named vectors (e.g. from q()) become a single-row table.
  #
  # gt renders to HTML, which is what you want in the RStudio viewer but not
  # at a terminal -- printing a gt object non-interactively dumps raw markup.
  # The parallel wrapper runs non-interactively, so fall back to a plain
  # data.frame print there.
  show_gt <- function(x, title, digits = 3, rowname_col = "parameter") {
    if (!is.data.frame(x)) {
      if (is.matrix(x) && !is.null(rownames(x))) {
        df <- as.data.frame(x)
        df <- cbind(setNames(data.frame(rownames(x)), rowname_col), df)
        rownames(df) <- NULL
        x <- df
      } else if (is.matrix(x)) {
        x <- as.data.frame(x)
      } else {
        x <- as.data.frame(as.list(x))
      }
    }
    num_cols <- names(x)[vapply(x, is.numeric, logical(1))]
    if (interactive()) {
      tbl <- gt::gt(x) |>
        gt::tab_header(title = title) |>
        gt::fmt_number(columns = num_cols, decimals = digits)
      print(tbl)
      return(invisible(tbl))
    }
    y <- x
    y[num_cols] <- lapply(y[num_cols], round, digits)
    cat("\n", title, "\n", strrep("-", nchar(title)), "\n", sep = "")
    print(y, row.names = FALSE)
    invisible(x)
  }

  report_note <- function(...) cat(paste0(...))

  # N_rec/Imm/N_surv/R_mean/S_mean are only defined for t >= 3 (recruitment
  # carries a two-year lag). Their t = 1,2 columns carry non-finite entries
  # that trip MCMCsummary's quantile call, so summarise only clean columns.
  .m <- as.matrix(samples)
  .ok <- colnames(.m)[!apply(.m, 2, function(x) any(!is.finite(x)))]
  print(MCMCsummary(samples[, .ok], round = 3))

  mat <- as.matrix(samples)
  q <- function(x) {
    c(
      mean = mean(x),
      lwr = quantile(x, .025, names = FALSE),
      upr = quantile(x, .975, names = FALSE)
    )
  }

  keyp <- intersect(
    c(
      "phi_pre",
      "phi_post",
      "pi_rec",
      "kappa",
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
  show_gt(
    round(gelman.diag(samples[, keyp], multivariate = FALSE)$psrf, 3),
    "Convergence (Rhat, key params; want < 1.05)"
  )

  pa <- t(apply(mat[, paste0("phi_ad[", 1:(T - 1), "]")], 2, q))
  show_gt(
    data.frame(
      transition = paste0(ipm_years[-T], "->", ipm_years[-1]),
      round(pa, 3)
    ),
    "Baseline adult survival phi_ad (untagged, female) by transition"
  )
  report_note(sprintf(
    "pooled pre-2014 (phi_pre): %.3f | pooled 2025->2026 (phi_post): %.3f\n",
    mean(mat[, "phi_pre"]),
    mean(mat[, "phi_post"])
  ))

  eff <- rbind(
    Sex_M_vs_F = q(mat[, "bSex_phi[2]"]),
    Sex_Unk_vs_F = q(mat[, "bSex_phi[3]"]),
    Tag_geo_vs_none = q(mat[, "bTag_phi[2]"]),
    Tag_GPS_vs_none = q(mat[, "bTag_phi[3]"])
  )
  show_gt(
    round(eff, 3),
    "Survival covariate effects (logit scale) and odds ratios",
    rowname_col = "effect"
  )
  report_note("Odds ratios (exp):\n")
  report_note(sprintf(
    "  geolocator vs untagged: %.2f (%.2f, %.2f)\n",
    exp(mean(mat[, "bTag_phi[2]"])),
    exp(quantile(mat[, "bTag_phi[2]"], .025)),
    exp(quantile(mat[, "bTag_phi[2]"], .975))
  ))
  report_note(sprintf(
    "  GPS vs untagged:        %.2f (%.2f, %.2f)\n",
    exp(mean(mat[, "bTag_phi[3]"])),
    exp(quantile(mat[, "bTag_phi[3]"], .025)),
    exp(quantile(mat[, "bTag_phi[3]"], .975))
  ))
  report_note(sprintf(
    "P(geolocator reduces survival): %.3f | P(GPS reduces survival): %.3f\n",
    mean(mat[, "bTag_phi[2]"] < 0),
    mean(mat[, "bTag_phi[3]"] < 0)
  ))

  report_note("\n--- Recapture: baseline + tag effects on p ---\n")
  report_note(sprintf(
    "baseline p (untagged): %.3f\n",
    plogis(mean(mat[, "lp_base"]))
  ))
  report_note(sprintf(
    "  tag effect on p (geo): OR %.2f | (GPS): OR %.2f\n",
    exp(mean(mat[, "bTag_p[2]"])),
    exp(mean(mat[, "bTag_p[3]"]))
  ))

  fr <- t(apply(mat[, paste0("f_rate[", 1:T, "]")], 2, q))
  show_gt(
    data.frame(year = ipm_years, round(fr, 3)),
    "Fecundity f_rate (chicks/pair) by year"
  )
  report_note(sprintf(
    "pooled 2026 (beyond count series): %.3f\n",
    mean(mat[, paste0("f_rate[", T + 1, "]")])
  ))

  lam <- t(apply(mat[, paste0("lambda[", 1:(T - 1), "]")], 2, q))
  show_gt(
    data.frame(
      transition = paste0(ipm_years[-T], "->", ipm_years[-1]),
      round(lam, 3)
    ),
    "Population growth lambda"
  )
  geolam <- exp(mean(log(rowMeans(mat[, paste0("lambda[", 1:(T - 1), "]")]))))
  report_note(sprintf(
    "Approx geometric-mean lambda: %.3f (%.1f%%/yr)\n",
    geolam,
    100 * (geolam - 1)
  ))

  rec <- rbind(
    pi_rec = q(mat[, "pi_rec"]),
    kappa = q(mat[, "kappa"]),
    psi1 = q(mat[, "psi1"]),
    psi2 = q(mat[, "psi2"])
  )
  show_gt(
    round(rec, 4),
    "Local recruitment (composite survival x philopatry)",
    digits = 4
  )
  # Priors are informative (1990s module); show how far the counts moved them.
  # If posterior ~= prior the count data carry little extra recruitment signal.
  pri <- rbind(
    pi_rec = c(
      mean = 2.268 / (2.268 + 13.610),
      lwr = qbeta(.025, 2.268, 13.610),
      upr = qbeta(.975, 2.268, 13.610)
    ),
    kappa = c(
      mean = 6.054 / (6.054 + 6.309),
      lwr = qbeta(.025, 6.054, 6.309),
      upr = qbeta(.975, 6.054, 6.309)
    )
  )
  cmp <- cbind(
    prior_mean = pri[, "mean"],
    post_mean = rec[rownames(pri), "mean"],
    prior_lwr = pri[, "lwr"],
    prior_upr = pri[, "upr"],
    post_lwr = rec[rownames(pri), "lwr"],
    post_upr = rec[rownames(pri), "upr"]
  )
  show_gt(
    round(cmp, 4),
    "Prior -> posterior movement (informative priors from 1990s module)",
    digits = 4
  )

  show_gt(round(q(mat[, "omega"]), 3), "Immigration omega")
  report_note(sprintf("P(omega > 0.01): %.3f\n", mean(mat[, "omega"] > 0.01)))
  # Recruitment and immigration both add birds; report their relative
  # contribution to annual gains so the trade-off is visible.
  Rm <- mat[, paste0("N_rec[", 3:T, "]"), drop = FALSE]
  Im <- mat[, paste0("Imm[", 3:T, "]"), drop = FALSE]
  frac <- rowSums(Rm) / (rowSums(Rm) + rowSums(Im))
  report_note(sprintf(
    "Share of annual gains from local recruitment: %.3f (%.3f, %.3f)\n",
    mean(frac),
    quantile(frac, .025),
    quantile(frac, .975)
  ))

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
  show_gt(
    group_survival(mat, 1),
    "Survival by sex x tag (transition 1; * GPS x Unknown has no data)",
    rowname_col = "group"
  )

  # ---- abundance trajectory plot ---------------------------------------------
  Nt <- t(apply(mat[, paste0("N_tot[", 1:T, "]")], 2, q))
  # png(file.path(plot_dir, "leye_ipm_abundance.png"), width=1600, height=1000, res=200)
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
  # dev.off()
  report_note("\n[figure] population trajectory drawn\n")

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

  # png(file.path(plot_dir, "leye_covariate_forest.png"), width = 1700, height = 1150, res = 200)
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
  # dev.off()
  report_note("[figure] covariate forest plot drawn\n\n")

  show_gt(
    round(phi_rows, 3),
    "Survival effects (odds-ratio scale)",
    rowname_col = "contrast"
  )
  show_gt(
    round(p_rows, 3),
    "Recapture effects (odds-ratio scale)",
    rowname_col = "contrast"
  )

  # =============================================================================
  # Figure: baseline adult annual survival from the LEYE IPM (covariate model)
  #
  # phi_ad[t] = baseline adult survival = untagged, FEMALE-reference birds.
  #             (This is the population rate driving the process model.)
  #
  # There is deliberately NO juvenile survival curve here. Juvenile survival is
  # no longer a model quantity: the old phi_juv = rho * phi_ad had rho weakly
  # identified (it traded off directly against immigration), and the 1990s chick
  # data showed recruitment is delayed and partly to age 2. Recruitment is now
  # psi1/psi2 -- composite survive-AND-return-locally probabilities, not annual
  # survival rates -- and is reported in the recruitment table above rather than
  # plotted on a survival axis, where it would invite a false comparison.
  #
  # Tagged/male survival differs by the fitted covariate effects; this figure
  # shows the BASELINE population rates. Use the forest plot for covariates.
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

  # transitions with direct CJS data (start year present in CJS intervals)
  has_cjs <- trans_yr %in% cjs_int_start

  # png(file.path(plot_dir, "leye_survival_adult.png"), width = 1700, height = 1050, res = 200)
  par(mar = c(4.5, 4.5, 2.5, 1))
  col_ad <- "#1f4e79"

  plot(
    NA,
    xlim = range(trans_yr) + c(-0.5, 0.5),
    ylim = c(0, 1),
    xlab = "Transition (survival from year t to t+1)",
    ylab = "Annual apparent survival",
    main = "LEYE baseline adult survival (untagged, female)",
    xaxt = "n"
  )
  axis(1, at = trans_yr, labels = paste0(trans_yr + 1), cex.axis = 1, las = 2)
  abline(h = seq(0, 1, .2), col = "grey92")

  arrows(
    trans_yr,
    ad[, "lwr"],
    trans_yr,
    ad[, "upr"],
    length = .025,
    angle = 90,
    code = 3,
    col = col_ad,
    lwd = 1.5
  )
  points(trans_yr, ad[, "med"], pch = 16, col = col_ad, cex = 1.2)
  lines(trans_yr, ad[, "med"], col = col_ad, lwd = 1, lty = 2)

  # mark transitions lacking direct CJS data (informed only via process model)
  # no_dat <- trans_yr[!has_cjs]
  # if (length(no_dat)) {
  #   mtext("*", side = 3, at = no_dat, line = -0.5, col = "grey40", cex = 1.5)
  # }

  mtext(
    # "Points = posterior median; bars = 95% CrI\n* = transition without direct CJS data (informed via process model).",
    "Points = posterior median; bars = 95% CrI",
    side = 3,
    line = -1.2,
    cex = 0.62,
    adj = 0.02
  )
  # dev.off()
  report_note("[figure] baseline adult survival drawn\n")

  show_gt(
    data.frame(transition = paste0(trans_yr, "-", trans_yr + 1), round(ad, 3)),
    "Adult survival (baseline)"
  )

  cat(
    "\n---------------------------------------------------------------\n",
    "Console summary complete. For the shareable report (figures\n",
    "embedded, interpretation and caveats included), run:\n\n",
    "    quarto render reports/leye_ipm_report.qmd\n",
    "---------------------------------------------------------------\n",
    sep = ""
  )

  invisible(list(mat = mat, phi_rows = phi_rows, p_rows = p_rows))
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

  # Was niter = nburnin = 10, which returns ZERO posterior draws and makes
  # every summary below fail. Use the parallel wrapper for production runs;
  # these settings make the serial path self-contained and actually runnable.
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

  ## ===========================================================================
  ## SECTION 5 -- SUMMARY
  ## (defined above; called here on this run's own 3-chain `samples`)
  ## ===========================================================================
  summary_out <- summarize_ipm(samples, ipm_years, T, cjs_int_start)
  mat <- summary_out$mat

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

  # ---------------------------------------------------------------------------
  # Retrospective correlations with population growth.
  #
  # NOTE: juvenile survival is no longer an annually-varying quantity, so the
  # old phi_juv-vs-lambda correlation has no analogue. Recruitment is now the
  # scalar pair psi1/psi2; what varies between years is the NUMBER of recruits
  # (N_rec) and immigrants (Imm) they generate. Those are the informative
  # annual series, so we correlate those instead.
  #
  # Index alignment: lambda[t] = N_tot[t+1]/N_tot[t], so growth over t->t+1 is
  # driven by gains realised IN year t+1, i.e. N_rec[t+1] and Imm[t+1].
  # N_rec/Imm exist for t = 3..T (recruitment carries a two-year lag), so the
  # gain series pairs with lambda[2..T-1].
  # ---------------------------------------------------------------------------
  lam_idx <- 2:(T - 1) # lambda indices pairing with gains in years 3..T
  corr <- matrix(NA, ncol = 4, nrow = n.draws)
  colnames(corr) <- c("phi_ad", "f_rate", "N_rec", "Imm")
  for (s in 1:n.draws) {
    corr[s, 1] <- cor(draws$phi_ad[s, 1:(T - 1)], draws$lambda[s, 1:(T - 1)])
    corr[s, 2] <- cor(draws$f_rate[s, 1:(T - 1)], draws$lambda[s, 1:(T - 1)])
    corr[s, 3] <- cor(draws$N_rec[s, 3:T], draws$lambda[s, lam_idx])
    corr[s, 4] <- cor(draws$Imm[s, 3:T], draws$lambda[s, lam_idx])
  }
  cri <- function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE)
  cat("\n--- Correlation of demographic series with lambda ---\n")
  # cor() is undefined for a draw whose series has no variance (e.g. Imm all
  # zero when omega ~ 0), so those draws drop out of the summary.
  print(round(
    cbind(mean = apply(corr, 2, mean, na.rm = TRUE), t(apply(corr, 2, cri))),
    3
  ))
  .nNA <- apply(corr, 2, function(x) sum(is.na(x)))
  if (any(.nNA > 0)) {
    cat(sprintf(
      "draws with an undefined (constant-series) correlation: %s of %d\n",
      paste(sprintf("%s=%d", names(.nNA), .nNA), collapse = ", "),
      n.draws
    ))
  }

  # ~~~~ Figure: lambda against the drivers ~~~~
  cri.lambda <- apply(draws$lambda, 2, cri)

  op <- par(mfrow = c(1, 3))
  layout(
    matrix(1:3, 1, 3, byrow = TRUE),
    widths = c(3.05, 3, 3),
    heights = rep(3.5, 3),
    TRUE
  )

  # panel 1: adult survival
  cri.rate <- apply(draws$phi_ad[, 1:(T - 1), drop = FALSE], 2, cri)
  plot(
    NA,
    ylim = range(cri.lambda),
    xlim = range(cri.rate),
    ylab = expression('Population growth rate (' * lambda * ')'),
    xlab = expression('Adult survival (' * phi[italic(a)] * ')'),
    axes = FALSE
  )
  axis(1)
  axis(2, las = 1)
  segments(
    colMeans(draws$phi_ad[, 1:(T - 1), drop = FALSE]),
    cri.lambda[1, ],
    colMeans(draws$phi_ad[, 1:(T - 1), drop = FALSE]),
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
  points(
    y = colMeans(draws$lambda),
    x = colMeans(draws$phi_ad[, 1:(T - 1), drop = FALSE]),
    pch = 16
  )

  # panel 2: productivity
  cri.rate <- apply(draws$f_rate[, 1:(T - 1), drop = FALSE], 2, cri)
  plot(
    NA,
    ylim = range(cri.lambda),
    xlim = range(cri.rate),
    ylab = NA,
    xlab = expression('Productivity (' * italic(f) * ')'),
    axes = FALSE
  )
  axis(1)
  axis(2, las = 1)
  segments(
    colMeans(draws$f_rate[, 1:(T - 1), drop = FALSE]),
    cri.lambda[1, ],
    colMeans(draws$f_rate[, 1:(T - 1), drop = FALSE]),
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
  points(
    y = colMeans(draws$lambda),
    x = colMeans(draws$f_rate[, 1:(T - 1), drop = FALSE]),
    pch = 16
  )

  # panel 3: recruits vs immigrants (the two sources of gains)
  rec_m <- colMeans(draws$N_rec[, 3:T, drop = FALSE])
  imm_m <- colMeans(draws$Imm[, 3:T, drop = FALSE])
  plot(
    NA,
    ylim = range(cri.lambda[, lam_idx]),
    xlim = range(c(rec_m, imm_m)),
    ylab = NA,
    xlab = "Annual gains (birds)",
    axes = FALSE
  )
  axis(1)
  axis(2, las = 1)
  points(
    y = colMeans(draws$lambda)[lam_idx],
    x = rec_m,
    pch = 16,
    col = "#1f4e79"
  )
  points(
    y = colMeans(draws$lambda)[lam_idx],
    x = imm_m,
    pch = 17,
    col = "#c0504d"
  )
  legend(
    "topleft",
    bty = "n",
    cex = 0.9,
    legend = c("local recruits", "immigrants"),
    col = c("#1f4e79", "#c0504d"),
    pch = c(16, 17)
  )
  par(op)

  # ---------------------------------------------------------------------------
  # LTRE -- DISABLED, NEEDS REWRITING FOR THE NEW STRUCTURE.
  #
  # The block that stood here was already non-functional before the
  # reparameterisation: it referenced draws$mean.f1, draws$mean.f2,
  # draws$mean.phij, draws$mean.phia, draws$mean.omega and wrote into
  # N.star.phij / N.star.om / N.star.f1 / N.star.f2, none of which are created
  # by this script (the arrays allocated just above it had different names).
  # It would have errored on the first iteration under the old model too.
  #
  # It also assumed a 3-stage projection matrix with two fecundity classes
  # (f1/f2) and recruitment at age 1, which no longer matches this model:
  # recruitment is now split across ages 1 and 2 via psi1/psi2 and there is a
  # single fecundity series. Rewriting it means deriving the transition matrix
  # for the delayed-recruitment structure, which is a separate task.
  #
  # To restore: build A = [[0, f*psi1/2, f*psi2/2], [phi_ad, 0, 0], ...] for the
  # age-structured population, then perturb psi1, psi2, phi_ad and omega in turn.
  # ---------------------------------------------------------------------------
} # end if(!SKIP_RUN)
