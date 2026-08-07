# =============================================================================
# LEYE transient LTRE (tLTRE) -- retrospective decomposition of population change
#
# Reads a cached IPM posterior; does NOT refit. Everything below is computed
# per posterior draw, so every contribution carries a credible interval rather
# than being a point estimate at posterior means.
#
# Two complementary questions, deliberately kept separate because they have
# different answers and are easy to conflate:
#
#   (A) LEVEL   -- what SUSTAINS the population?
#                  Decomposes mean lambda into survival, local recruitment and
#                  immigration. Exact and additive.
#
#   (B) VARIANCE (the tLTRE proper) -- what DRIVES year-to-year variation?
#                  Decomposes Var(lambda) across time into contributions from
#                  the demographic drivers, following Koons et al. (2016, 2017):
#                    contribution_i = sum_j Cov(theta_i, theta_j) * s_i * s_j
#                  with sensitivities s evaluated at each draw's temporal means.
#
# ---------------------------------------------------------------------------
# THE DECOMPOSITION
#
# From the process model (scripts/LEYE_IPM.R, section 2b), with s indexing the
# transition lambda[s] = N[s+1]/N[s]:
#
#   N[s+1] = N_surv[s+1] + N_rec[s+1] + Imm[s+1]
#
# Dividing through by N[s] gives an EXACT realised identity:
#
#   lambda[s] = surv_pc[s] + rec_pc[s] + imm_pc[s]                        (A)
#
# and taking expectations gives the driver-level form:
#
#   E[lambda[s]] = phi[s] + omega + (psi1/2)*f[s] + (psi2/2)*rho[s]*f[s-1]  (B)
#
# where rho[s] = N[s-1]/N[s] carries the delayed-recruitment population
# structure -- the transient term that an asymptotic LTRE would discard.
#
# ---------------------------------------------------------------------------
# TWO THINGS TO KNOW BEFORE READING THE OUTPUT
#
# 1. psi1, psi2 and omega are SCALARS in the current model -- they have no
#    year index. Their temporal variance is exactly zero, so their tLTRE
#    contribution is zero BY CONSTRUCTION, not by evidence. This script
#    reports them explicitly as structural zeros rather than silently
#    omitting them. If you want to ask whether immigration drives variation
#    in growth, omega must first be made time-varying in the IPM.
#
# 2. This population is small (N ~ 25), so Poisson demographic stochasticity
#    is a first-order term, not a rounding error. Identity (A) is exact but
#    includes that noise; (B) is the smooth expectation. The script partitions
#    Var(lambda) into a driver-driven part and a demographic-noise part, and
#    reports both. Attributing demographic noise to a "driver" would overstate
#    how much biology is being explained.
# =============================================================================

library(coda)
set.seed(20260805)

## ===========================================================================
## SECTION 1 -- LOAD
## ===========================================================================
tltre_load <- function(samples_path = NULL) {
  # Canonical filename, NOT a newest-match glob. A glob on LEYE_IPM_samples*.rds
  # also matches companion and archive fits (vague-priors, time-varying-omega,
  # dated caches); taking whichever is newest silently decomposes whichever
  # model happened to be refitted last.
  if (is.null(samples_path)) {
    samples_path <- "LEYE_IPM_samples.rds"
    if (!file.exists(samples_path)) {
      stop(
        "No cached posterior at '", samples_path, "'. Fit the model first with ",
        "scripts/LEYE_IPM_run_parallel.R (or scripts/LEYE_IPM.R), or pass ",
        "samples_path explicitly."
      )
    }
  }
  s <- readRDS(samples_path)
  # The parallel wrapper's lapply() used to strip the mcmc.list class, and
  # as.matrix() on a plain list silently returns a 3 x 1 matrix. Guard.
  if (!inherits(s, "mcmc.list")) s <- coda::as.mcmc.list(s)
  m <- as.matrix(s)
  if (nrow(m) < 100) {
    stop("Only ", nrow(m), " draws recovered -- the cache is probably corrupt.")
  }
  attr(m, "path") <- samples_path
  m
}

# Year axis, taken from the count data exactly as LEYE_IPM.R defines it, so
# this script stays consistent without having to source the whole model.
tltre_years <- function() {
  dc <- read.csv("./Arin LEYE R code/data_2025.csv")
  sort(unique(dc$year))
}

## ===========================================================================
## SECTION 2 -- PER-DRAW ANNUAL SERIES
## ===========================================================================
# Returns a list of (draws x transitions) matrices. Transitions are the s for
# which the identity is defined: s+1 must be >= 3 (recruitment's two-year lag)
# and s+1 <= T, so s runs 2..T-1.
tltre_series <- function(m, Tt) {
  g <- function(stem, i) m[, paste0(stem, "[", i, "]")]
  S <- 2:(Tt - 1)
  grab <- function(f) vapply(S, f, numeric(nrow(m)))

  N_s <- grab(function(i) g("N_ad", i))

  # omega may be a scalar (older fits) or a year-indexed vector. Carry it as a
  # matrix either way: if it is constant its temporal variance is zero and the
  # tLTRE returns a zero contribution automatically, with no special-casing.
  # NOTE the index. omega[t] is the rate for the year immigrants ARRIVE, so the
  # rate acting on transition s is omega[s+1], not omega[s].
  omega_tv <- any(grepl("^omega\\[", colnames(m)))
  omega_m <- if (omega_tv) {
    grab(function(i) g("omega", i + 1))
  } else {
    matrix(m[, "omega"], nrow(m), length(S))
  }

  out <- list(
    S = S,
    omega_time_varying = omega_tv,
    lambda = grab(function(i) g("lambda", i)),
    # (A) exact realised per-capita pathways
    surv_pc = grab(function(i) g("N_surv", i + 1)) / N_s,
    rec_pc = grab(function(i) g("N_rec", i + 1)) / N_s,
    imm_pc = grab(function(i) g("Imm", i + 1)) / N_s,
    # (B) drivers
    phi = grab(function(i) g("phi_ad", i)),
    f_t = grab(function(i) g("f_rate", i)),
    f_lag = grab(function(i) g("f_rate", i - 1)),
    rho = grab(function(i) g("N_ad", i - 1)) / N_s,
    omega = omega_m,
    # scalars
    psi1 = m[, "psi1"],
    psi2 = m[, "psi2"]
  )
  out$lambda_exp <- out$phi + out$omega +
    (out$psi1 / 2) * out$f_t +
    (out$psi2 / 2) * out$rho * out$f_lag
  out
}

# Fail loudly if the process model changes underneath this script. The exact
# identity (A) is the foundation of everything below; if it stops holding, the
# decomposition is silently wrong rather than obviously broken.
tltre_check <- function(x, tol = 1e-8) {
  err <- max(abs(x$lambda - (x$surv_pc + x$rec_pc + x$imm_pc)))
  if (err > tol) {
    stop(
      "Realised decomposition no longer exact (max error ", format(err),
      "). The process model in scripts/LEYE_IPM.R has changed -- rederive ",
      "the decomposition in this script's header before trusting any output."
    )
  }
  invisible(err)
}

## ===========================================================================
## SECTION 3 -- (A) LEVEL: what sustains the population?
## ===========================================================================
# Mean lambda splits additively across the three pathways. This is the
# "where do birds come from" question, and is NOT what the variance
# decomposition answers.
tltre_level <- function(x) {
  comp <- list(
    Survival = rowMeans(x$surv_pc),
    Recruitment = rowMeans(x$rec_pc),
    Immigration = rowMeans(x$imm_pc)
  )
  lam <- rowMeans(x$lambda)
  out <- t(vapply(comp, function(v) {
    c(
      contribution = mean(v),
      lwr = quantile(v, .025, names = FALSE),
      upr = quantile(v, .975, names = FALSE),
      pct_of_lambda = 100 * mean(v / lam)
    )
  }, numeric(4)))
  attr(out, "lambda") <- c(
    mean = mean(lam),
    lwr = quantile(lam, .025, names = FALSE),
    upr = quantile(lam, .975, names = FALSE)
  )
  out
}

## ===========================================================================
## SECTION 4 -- (B) VARIANCE: the tLTRE proper
## ===========================================================================
# Koons-style contributions, computed within each posterior draw:
#   contribution_i = sum_j Cov(theta_i, theta_j) * s_i * s_j
# Sensitivities are the partial derivatives of (B) w.r.t. each driver,
# evaluated at that draw's temporal means.
#
# Drivers carried: phi, f_t, f_lag, rho. psi1/psi2/omega are scalars and are
# reported separately as structural zeros.
tltre_variance <- function(x) {
  nd <- nrow(x$lambda)
  drivers <- c("phi", "omega", "f_t", "f_lag", "rho")
  contrib <- matrix(NA_real_, nd, length(drivers), dimnames = list(NULL, drivers))
  var_exp <- numeric(nd)
  var_tot <- numeric(nd)
  var_dem <- numeric(nd)

  cov_cross <- numeric(nd)
  lin_resid <- numeric(nd)

  for (d in seq_len(nd)) {
    Th <- cbind(
      phi = x$phi[d, ],
      omega = x$omega[d, ],
      f_t = x$f_t[d, ],
      f_lag = x$f_lag[d, ],
      rho = x$rho[d, ]
    )
    mu <- colMeans(Th)
    # partial derivatives of (B), at this draw's temporal means
    sens <- c(
      phi = 1,
      omega = 1,
      f_t = x$psi1[d] / 2,
      f_lag = (x$psi2[d] / 2) * mu[["rho"]],
      rho = (x$psi2[d] / 2) * mu[["f_lag"]]
    )
    C <- stats::cov(Th)
    contrib[d, ] <- as.vector(C %*% sens) * sens
    var_exp[d] <- stats::var(x$lambda_exp[d, ])
    var_tot[d] <- stats::var(x$lambda[d, ])
    resid <- x$lambda[d, ] - x$lambda_exp[d, ]
    var_dem[d] <- stats::var(resid)
    # lambda_exp and the demographic residual are NOT independent: rho[s]
    # depends on realised N, which carries the same Poisson noise. So the
    # partition needs its cross term or the shares can exceed 100%.
    cov_cross[d] <- 2 * stats::cov(x$lambda_exp[d, ], resid)
    # tLTRE is a first-order approximation and (B) contains a product term
    # (rho * f_lag), so the contributions need not sum exactly to var_exp.
    lin_resid[d] <- var_exp[d] - sum(contrib[d, ])
  }

  q3 <- function(v) {
    c(
      mean = mean(v),
      lwr = quantile(v, .025, names = FALSE),
      upr = quantile(v, .975, names = FALSE)
    )
  }
  tab <- rbind(
    t(apply(contrib, 2, q3)),
    psi1 = c(0, 0, 0),
    psi2 = c(0, 0, 0)
  )
  # share of the DRIVER-EXPLAINED variance (i.e. of var_exp, not var_tot)
  share <- 100 * colMeans(contrib / pmax(var_exp, .Machine$double.eps))
  tab <- cbind(tab, pct_of_explained = c(share, 0, 0))

  # Exact partition of Var(lambda). These three sum to the total by
  # construction, so their shares sum to 100%.
  pct <- function(v) q3(100 * v / var_tot)

  list(
    contributions = tab,
    variance = rbind(
      total = q3(var_tot),
      driver_explained = q3(var_exp),
      demographic = q3(var_dem),
      cross_term = q3(cov_cross)
    ),
    partition_pct = rbind(
      driver_explained = pct(var_exp),
      demographic = pct(var_dem),
      cross_term = pct(cov_cross)
    ),
    partition_check = q3(var_tot - (var_exp + var_dem + cov_cross)),
    linearisation_resid = q3(lin_resid)
  )
}

## ===========================================================================
## SECTION 5 -- RUN
## ===========================================================================
if (!exists("SKIP_TLTRE_RUN")) SKIP_TLTRE_RUN <- FALSE
if (!SKIP_TLTRE_RUN) {
  mat <- tltre_load()
  ipm_years <- tltre_years()
  Tt <- length(ipm_years)
  x <- tltre_series(mat, Tt)
  err <- tltre_check(x)

  cat("=====================================================================\n")
  cat("LEYE transient LTRE\n")
  cat("=====================================================================\n")
  cat(sprintf("Posterior : %s\n", basename(attr(mat, "path"))))
  cat(sprintf("Draws     : %s\n", format(nrow(mat), big.mark = ",")))
  cat(sprintf(
    "Transitions: %d (%d-%d to %d-%d)\n",
    length(x$S),
    ipm_years[min(x$S)], ipm_years[min(x$S) + 1],
    ipm_years[max(x$S)], ipm_years[max(x$S) + 1]
  ))
  cat(sprintf("Identity check: max |lambda - (surv+rec+imm)| = %.2e\n\n", err))

  lev <- tltre_level(x)
  lam <- attr(lev, "lambda")
  cat("---------------------------------------------------------------------\n")
  cat("(A) LEVEL -- what sustains the population?\n")
  cat("    Mean lambda splits additively across pathways.\n")
  cat("---------------------------------------------------------------------\n")
  cat(sprintf(
    "mean lambda = %.4f (%.4f, %.4f)\n\n", lam[1], lam[2], lam[3]
  ))
  print(round(lev, 4))

  vr <- tltre_variance(x)
  cat("\n---------------------------------------------------------------------\n")
  cat("(B) VARIANCE -- what drives year-to-year variation?\n")
  cat("---------------------------------------------------------------------\n")
  cat("\nVariance of lambda across transitions:\n")
  print(round(vr$variance, 6))
  cat(sprintf(
    "  (partition closes to %.2e -- total = explained + demographic + cross)\n",
    abs(vr$partition_check[1])
  ))
  cat("\nShare of Var(lambda), summing to 100%:\n")
  print(round(vr$partition_pct, 1))
  cat("\nThe population is small (N ~ 25), so Poisson noise is a real term, not\n")
  cat("an artefact. The cross term is non-zero because rho depends on realised\n")
  cat("N and so carries the same noise -- the smooth and stochastic parts of\n")
  cat("lambda are not independent. Because that term can be negative, an\n")
  cat("individual share may exceed 100% even though the three always sum to it.\n")

  cat("\ntLTRE contributions to the DRIVER-EXPLAINED variance:\n")
  print(round(vr$contributions, 6))
  cat(sprintf(
    "  (linearisation residual %.2e; tLTRE is first-order and (B) has a\n   rho*f_lag product term, so contributions need not sum exactly)\n",
    vr$linearisation_resid[1]
  ))
  cat("\nNOTE: psi1 and psi2 are scalars in the IPM, so their zero contribution\n")
  cat("is structural -- they cannot vary between years by construction.\n")
  if (isTRUE(x$omega_time_varying)) {
    cat("\nomega IS time-varying in this fit, so its contribution is estimated.\n")
    cat("Read it against the demographic share above: with ~3 immigrants a\n")
    cat("year, a year effect can absorb Poisson noise and re-present it as\n")
    cat("environmental variation. If the driver-explained share rose while the\n")
    cat("demographic share fell by a matching amount relative to the\n")
    cat("constant-omega fit, that is relabelling, not new signal. Check\n")
    cat("sd_omega against its prior before interpreting.\n")
  } else {
    cat("omega is also a scalar here, so it too contributes zero by\n")
    cat("construction. Make it time-varying to ask whether immigration drives\n")
    cat("variation in growth.\n")
  }

  res <- list(
    level = lev,
    variance = vr,
    series = x,
    years = ipm_years,
    samples_path = attr(mat, "path"),
    generated = Sys.time()
  )
  saveRDS(res, "LEYE_tLTRE_results.rds")
  cat("\nWrote: LEYE_tLTRE_results.rds\n")
}
