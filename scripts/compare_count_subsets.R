# =============================================================================
# Count-stream sensitivity: does the survey panel shape the trend?
#
# Run from the PROJECT ROOT, after the companion fits exist:
#     Rscript scripts/LEYE_IPM_sensitivity.R full
#     Rscript scripts/LEYE_IPM_sensitivity.R balanced9
#     Rscript scripts/compare_count_subsets.R
#
# THE QUESTION
# The IPM's three data streams did not cover the same ground. Every marked bird
# and every monitored nest comes from the Anchorage bowl, while the counts came
# from 32 wetlands spanning 67 km -- and five of those, 37-67 km out, supplied
# 22.6% of all birds counted. On top of that the panel is unstable: it grew
# from 13 wetlands in 2014 to 32 by 2021, only 9 were surveyed throughout, and
# the wetlands added late are the rich ones (correlation between a site's mean
# count and its mean survey year, +0.50).
#
# Since 2026-11-04 the MAIN model drops those five wetlands, so all three
# streams describe the Anchorage bowl. This script is the record of why, and
# the standing check that the choice is not doing something unintended:
#
#   main       the Anchorage bowl, 780 records -- LEYE_IPM_samples.rds
#   full       all 32 wetlands, 935 records: what the wider area does
#   balanced9  the 9 wetlands surveyed in all 13 years, 333 records: whether
#              the changing panel manufactures trend on its own
#
# All three share priors, chain length, seeds and monitors, so the only
# difference is which wetlands enter the count likelihood.
#
# WHAT TO READ, AND WHAT TO IGNORE
# Read lambda, P(lambda<1), and the shape of N_tot. Ignore logB0 and the
# absolute level of N_tot: logB0 is weakly identified against the population
# level in every fit (cor about -0.85 with mean log N_ad, ESS ~100 of 21000),
# so it differs between fits for reasons that have nothing to do with the
# subset.
#
# OUTPUT
#   console table, and output/plots/count_subset_comparison.png
# =============================================================================

suppressMessages({
  library(coda)
})

FITS <- c(
  main = "LEYE_IPM_samples.rds",
  full = "LEYE_IPM_samples_full.rds",
  balanced9 = "LEYE_IPM_samples_balanced9.rds"
)
LABEL <- c(
  main = "27, Anchorage bowl",
  full = "all 32 wetlands",
  balanced9 = "9, surveyed all years"
)

missing <- FITS[!file.exists(FITS)]
if (length(missing)) {
  stop(
    "Missing fit(s):\n  ", paste(missing, collapse = "\n  "),
    "\nRun scripts/LEYE_IPM_sensitivity.R with the matching setting name(s)."
  )
}
sv <- lapply(FITS, readRDS)

## ---------------------------------------------------------------------------
## Summaries
## ---------------------------------------------------------------------------
# GEOMETRIC mean across years, within each draw -- the growth rate, and the
# summary the report quotes. rowMeans() of lambda is the arithmetic mean and is
# not a growth rate; see the headline chunk of reports/leye_ipm_report.qmd.
geolam <- function(m) {
  lc <- grep("lambda[", colnames(m), fixed = TRUE, value = TRUE)
  exp(rowMeans(log(m[, lc, drop = FALSE])))
}
ntot <- function(m) m[, grep("N_tot[", colnames(m), fixed = TRUE, value = TRUE), drop = FALSE]

MATS <- lapply(sv, as.matrix)
ci <- function(x) c(mean(x), quantile(x, .025), quantile(x, .975))

cat("\n=====================================================================\n")
cat("COUNT-STREAM SENSITIVITY\n")
cat("=====================================================================\n")
cat(sprintf(
  "%-11s %-22s %-24s %7s %9s %9s\n",
  "setting", "count stream", "geometric-mean lambda", "P(<1)", "%/yr", "N_T/N_1"
))
for (k in names(FITS)) {
  m <- MATS[[k]]
  g <- geolam(m)
  n <- ntot(m)
  tot <- n[, ncol(n)] / n[, 1]
  cat(sprintf(
    "%-11s %-22s %.3f (%.3f-%.3f)%s %7.3f %8.1f%% %8.3f\n",
    k, LABEL[k], mean(g), quantile(g, .025), quantile(g, .975),
    "    ", mean(g < 1), 100 * (mean(g) - 1), mean(tot)
  ))
}

# Demographic parameters. These share their data across all three fits, so a
# shift here means the count stream is pulling on them through the shared
# latent series rather than that the data changed.
cat("\n--- demographic parameters (same data in every fit) ---\n")
PARS <- c("pi_rec", "kappa", "omega", "sd_site", "sd_od")
cat(sprintf("%-11s %s\n", "setting", paste(sprintf("%-20s", PARS), collapse = "")))
for (k in names(FITS)) {
  m <- MATS[[k]]
  v <- vapply(PARS, function(p) {
    if (!p %in% colnames(m)) return("--")
    x <- ci(m[, p])
    sprintf("%.3f (%.3f-%.3f)", x[1], x[2], x[3])
  }, character(1))
  cat(sprintf("%-11s %s\n", k, paste(sprintf("%-20s", v), collapse = "")))
}

# Baseline adult survival, averaged over years: the quantity the mark-resight
# data pin directly, and therefore the one that should move LEAST.
cat("\n--- mean baseline adult survival (untagged female) ---\n")
for (k in names(FITS)) {
  m <- MATS[[k]]
  ph <- m[, grep("phi_ad[", colnames(m), fixed = TRUE, value = TRUE), drop = FALSE]
  x <- ci(rowMeans(ph))
  cat(sprintf("  %-11s %.3f (%.3f-%.3f)\n", k, x[1], x[2], x[3]))
}

cat("\n--- convergence (max Rhat over lambda, N_tot and the parameters above) ---\n")
for (k in names(FITS)) {
  p <- grep("^(lambda|N_tot)\\[|^(pi_rec|kappa|omega|sd_site|sd_od|logB0)$",
            colnames(sv[[k]][[1]]), value = TRUE)
  r <- gelman.diag(sv[[k]][, p], multivariate = FALSE)$psrf[, 1]
  cat(sprintf("  %-11s max Rhat %.3f (%s)%s\n", k, max(r), names(which.max(r)),
              if (max(r) > 1.05) "   <- CHECK" else ""))
}

## ---------------------------------------------------------------------------
## Does site_re track observed richness, or is it shrunk?
## ---------------------------------------------------------------------------
# The reason for monitoring site_re at all. If a rich wetland's fitted level is
# pulled toward the mean, the surplus has nowhere to go but N_tot for the years
# that wetland was surveyed -- and the rich wetlands here were surveyed late,
# so the leak would flatten the trend and understate the decline.
#
# A slope of 1 on log(1 + mean count) means no shrinkage at all. Well under 1,
# concentrated in the wetlands with fewest survey-years, is the failure mode.
dc <- read.csv("Arin LEYE R code/data_2026.csv", fileEncoding = "UTF-8-BOM")
pres <- table(dc$wetland, dc$year) > 0
site_tab <- function(k, keep) {
  m <- MATS[[k]]
  sc <- grep("site_re[", colnames(m), fixed = TRUE, value = TRUE)
  if (!length(sc)) return(NULL)
  ids <- sort(unique(dc$wetland[dc$wetland %in% keep]))
  stopifnot(length(ids) == length(sc)) # site_re[j] is the j-th wetland, sorted
  d <- dc[dc$wetland %in% keep, ]
  data.frame(
    id = ids,
    site_re = colMeans(m[, sc]),
    raw = as.numeric(tapply(d$count, d$wetland, mean)[as.character(ids)]),
    yrs = as.integer(rowSums(pres)[as.character(ids)]),
    meanyr = as.numeric(tapply(d$year, d$wetland, mean)[as.character(ids)])
  )
}
KEEP <- list(
  main = setdiff(sort(unique(dc$wetland)), c(3, 10, 18, 22, 23)),
  full = sort(unique(dc$wetland)),
  balanced9 = as.integer(rownames(pres)[rowSums(pres) == ncol(pres)])
)
st <- lapply(names(FITS), function(k) site_tab(k, KEEP[[k]]))
names(st) <- names(FITS)

if (!is.null(st$full)) {
  cat("\n--- site_re against observed richness ---\n")
  for (k in names(FITS)) {
    d <- st[[k]]
    if (is.null(d)) {
      cat(sprintf("  %-11s site_re not monitored in this fit\n", k))
      next
    }
    fit <- lm(site_re ~ log1p(raw), d)
    cat(sprintf(
      "  %-11s slope %.2f   cor %.3f   cor(site_re, mean survey year) %+.3f\n",
      k, coef(fit)[2], cor(d$site_re, log1p(d$raw)), cor(d$site_re, d$meanyr)
    ))
  }
  d <- st$full
  d$resid <- resid(lm(site_re ~ log1p(raw), d))
  d <- d[order(d$resid), ]
  cat("\n  Wetlands whose fitted level sits furthest BELOW their richness\n")
  cat("  (the direction that leaks abundance into N_tot):\n")
  show <- head(d, 5)
  for (i in seq_len(nrow(show))) {
    cat(sprintf(
      "    id %2d  raw %5.2f/survey  %2d yrs  mean yr %.1f  site_re %+.2f  resid %+.2f\n",
      show$id[i], show$raw[i], show$yrs[i], show$meanyr[i],
      show$site_re[i], show$resid[i]
    ))
  }
}

## ---------------------------------------------------------------------------
## Figure
## ---------------------------------------------------------------------------
VER <- "#D55E00"; BLU <- "#0072B2"; GRY <- "#8A96A0"
DECK_BG <- "#FAF8F4"; INK <- "#1C2B33"; MID <- "#5A6872"
COLS <- c(main = INK, full = BLU, balanced9 = VER)
fade <- function(h, a) {
  v <- col2rgb(h)
  rgb(v[1], v[2], v[3], alpha = a * 255, maxColorValue = 255)
}

yrs <- sort(unique(dc$year))
dir.create("output/plots", showWarnings = FALSE, recursive = TRUE)
png("output/plots/count_subset_comparison.png",
    width = 1300, height = 560, res = 110, type = "cairo", bg = DECK_BG)
par(mfrow = c(1, 2), bg = DECK_BG, family = "sans", col.axis = MID,
    col.lab = INK, fg = MID, mar = c(4.4, 4.8, 3.4, 1.2))

# N_tot is a RELATIVE index, so the panels are normalised to each fit's own
# first year. Plotting the raw level would compare three arbitrary scalings of
# logB0 and show nothing about trend.
plot(NA, xlim = range(yrs), ylim = c(0, 1.5), xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "Breeding adults, relative to 2014")
abline(h = seq(0, 1.5, .25), col = "#DCD6CB", lwd = 1.3)
for (k in names(FITS)) {
  n <- ntot(MATS[[k]])
  rel <- n / n[, 1]
  q <- apply(rel, 2, quantile, c(.025, .5, .975))
  polygon(c(yrs, rev(yrs)), c(q[1, ], rev(q[3, ])), col = fade(COLS[k], .13), border = NA)
  lines(yrs, q[2, ], col = COLS[k], lwd = 3)
}
axis(1, yrs[seq(1, length(yrs), 2)], lwd = 0, lwd.ticks = 1)
axis(2, seq(0, 1.5, .25), las = 1, lwd = 0, lwd.ticks = 1)
mtext("Fitted trajectory", side = 3, line = 1.4, adj = 0, cex = 1.15, font = 2, col = INK)
mtext("median and 95% CRI, each normalised to its own 2014",
      side = 3, line = 0.2, adj = 0, cex = .92, col = MID)
legend("bottomleft", legend = sprintf("%s (%s)", names(FITS), LABEL[names(FITS)]),
       col = COLS[names(FITS)], lwd = 3, bty = "n", cex = .92, text.col = INK)

# Growth rate: the number the manuscript quotes.
g <- lapply(MATS, geolam)
xr <- range(unlist(lapply(g, quantile, c(.005, .995))))
plot(NA, xlim = xr, ylim = c(0.5, length(g) + 0.5), yaxt = "n", xaxt = "n",
     bty = "n", xlab = "Geometric-mean annual lambda", ylab = "")
abline(v = 1, lty = 2, lwd = 2, col = INK)
for (i in seq_along(g)) {
  k <- names(g)[i]
  q <- quantile(g[[k]], c(.025, .975))
  y <- length(g) - i + 1
  arrows(q[1], y, q[2], y, code = 3, angle = 90, length = .07, lwd = 4, col = COLS[k])
  points(mean(g[[k]]), y, pch = 19, cex = 2.1, col = COLS[k])
  text(mean(g[[k]]), y + .22, sprintf("%.3f   P(<1) = %.3f", mean(g[[k]]), mean(g[[k]] < 1)),
       cex = .92, font = 2, col = COLS[k])
}
axis(1, lwd = 0, lwd.ticks = 1)
axis(2, at = length(g):1, labels = names(g), las = 1, lwd = 0, cex.axis = .95,
     col.axis = INK)
mtext("Growth rate", side = 3, line = 1.4, adj = 0, cex = 1.15, font = 2, col = INK)
mtext("dashed line is a stable population", side = 3, line = 0.2, adj = 0,
      cex = .92, col = MID)
invisible(dev.off())

cat("\nWrote: output/plots/count_subset_comparison.png\n")
