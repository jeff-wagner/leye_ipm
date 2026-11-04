# Recruitment / local return rate: 1990s vs modern era.
#
# Run from the PROJECT ROOT (all paths are relative, as elsewhere in scripts/):
#     Rscript scripts/fig_era_recruitment.R
#
# Reads LEYE_recruitment_two_era.rds and writes six variants of one figure
# into manuscript/figures/ -- all from a single body(), so a deck version and
# its talk export cannot drift apart:
#
#   era_recruitment.png            cohort points, deck cream background
#   era_recruitment_simple.png     era estimates only, deck cream background
#   era_recruitment_talk.png/.pdf  cohort points, TRANSPARENT, for a talk
#   era_recruitment_simple_talk.*  era estimates only, TRANSPARENT, for a talk
#
# The transparent variants drop onto any LIGHT slide template. Text is dark
# ink, so they need the palette inverted for a dark background.
#
# The report draws this same contrast live from the module posterior rather
# than embedding these files, so the two cannot disagree.
#
# SENSITIVITY OVERLAY
# When LEYE_recruitment_two_era_drop1997.rds is present (written by
# scripts/LEYE_recruitment_drop_cohort.R), both panels also show the fit with
# the 1997 cohort removed: hollow markers on the left, a dashed posterior on
# the right. That cohort returned 7 of 11 chicks, an order of magnitude above
# every other cohort in either era, and it sets the 1990s baseline the whole
# contrast is measured against -- so the figure should not imply the decline
# rests on more than it does.
#
# Set SHOW_DROP97 <- FALSE to restore the plain two-era figure exactly.
#
m <- do.call(rbind, readRDS("LEYE_recruitment_two_era.rds"))
pic <- grep("pi_c[", colnames(m), fixed = TRUE, value = TRUE)
N90 <- 5                      # 1995-1999 cohorts; modern cohorts follow
c90 <- colMeans(m[, pic[1:N90]]); cM <- colMeans(m[, pic[(N90 + 1):length(pic)]])
coh90 <- 1995:1999

DROP_PATH <- "LEYE_recruitment_two_era_drop1997.rds"
SHOW_DROP97 <- file.exists(DROP_PATH)

VER <- "#D55E00"; BLU <- "#0072B2"; GRY <- "#8A96A0"
DECK_BG <- "#FAF8F4"; INK <- "#1C2B33"; MID <- "#5A6872"
fade <- function(h, a) { v <- col2rgb(h); rgb(v[1], v[2], v[3], alpha = a * 255, maxColorValue = 255) }
ci <- function(x) c(mean(x), quantile(x, .025), quantile(x, .975))

raw <- function(f, occ) {
  d <- read.csv(f, check.names = FALSE)
  eh <- as.matrix(d[, occ]); eh[is.na(eh)] <- 0
  c(ret = sum(apply(eh, 1, function(x) sum(x == 1) > 1)), n = nrow(d))
}
r90 <- raw("data/LEYE_juv_EH_1995_2000.csv", as.character(1995:2000))
rM  <- raw("data/LEYE_chick_EH_2018_2026.csv", as.character(2018:2026))
e90 <- ci(m[, "pi_90"]); eM <- ci(m[, "pi_M"])
rr <- m[, "era_ratio"]

if (SHOW_DROP97) {
  md <- do.call(rbind, readRDS(DROP_PATH))
  e90d <- ci(md[, "pi_90"]); eMd <- ci(md[, "pi_M"])
  rrd <- md[, "era_ratio"]
}

body <- function(show_cohorts, bg) {
  layout(matrix(1:2, 1), widths = c(1.15, 1))
  par(bg = bg, col.axis = MID, col.lab = INK, fg = MID, family = "sans",
      cex.axis = 1.2, cex.lab = 1.35, lend = 1)
  hollow <- if (identical(bg, "transparent")) NA else bg

  par(mar = c(5.0, 6.2, 3.4, 1.0))
  ytop <- if (show_cohorts) 0.65 else 0.38
  step <- if (show_cohorts) 0.1 else 0.05
  plot(NA, xlim = c(0.45, 2.55), ylim = c(0, ytop), xaxt = "n", yaxt = "n",
       bty = "n", xlab = "", ylab = "Probability of surviving and returning to breed")
  abline(h = seq(0, ytop, step), col = "#DCD6CB", lwd = 1.4)
  if (show_cohorts) {
    set.seed(1)
    x90 <- jitter(rep(1, length(c90)), amount = .085)
    points(x90, c90, pch = 21, bg = fade(VER, .40), col = "white", lwd = 1.4, cex = 1.5)
    points(jitter(rep(2, length(cM)), amount = .085), cM,
           pch = 21, bg = fade(BLU, .40), col = "white", lwd = 1.4, cex = 1.5)
    hi <- which.max(c90)
    if (SHOW_DROP97) {
      points(x90[hi], c90[hi], pch = 21, bg = NA, col = VER, lwd = 2.6, cex = 2.2)
    }
    text(x90[hi], c90[hi], sprintf("%d cohort", coh90[hi]), pos = 4,
         offset = 0.7, cex = 1.1, col = VER)
  }
  # sensitivity fit sits to the LEFT of each era, clear of the value labels
  if (SHOW_DROP97) {
    xo <- c(1, 2) - 0.22
    arrows(xo, c(e90d[2], eMd[2]), xo, c(e90d[3], eMd[3]),
           code = 3, angle = 90, length = .07, lwd = 2.4, col = c(VER, BLU))
    points(xo, c(e90d[1], eMd[1]), pch = 21, bg = hollow,
           col = c(VER, BLU), lwd = 2.8, cex = 1.9)
  }
  arrows(c(1, 2), c(e90[2], eM[2]), c(1, 2), c(e90[3], eM[3]),
         code = 3, angle = 90, length = .10, lwd = 4.5, col = c(VER, BLU))
  points(c(1, 2), c(e90[1], eM[1]), pch = 19, cex = 2.6, col = c(VER, BLU))
  text(c(1, 2), c(e90[1], eM[1]), sprintf("%.3f", c(e90[1], eM[1])),
       pos = 4, offset = 1.2, cex = 1.35, font = 2, col = c(VER, BLU))
  axis(1, 1:2, c("1995\u20131999", "2018\u20132025"), lwd = 0, cex.axis = 1.35, col.axis = INK)
  axis(2, seq(0, ytop, step), las = 1, lwd = 0, lwd.ticks = 1)
  mtext("Modeled recruitment, with 95% CRI", side = 3, line = 1.4, adj = 0,
        cex = 1.35, font = 2, col = INK)
  sub <- sprintf("raw returns %d/%d vs %d/%d", r90["ret"], r90["n"], rM["ret"], rM["n"])
  if (show_cohorts) sub <- paste("small points = individual cohorts   |  ", sub)
  if (SHOW_DROP97) sub <- paste(sub, "  |   hollow = 1997 cohort excluded")
  mtext(sub, side = 3, line = 0.1, adj = 0, cex = 1.05, col = MID)

  par(mar = c(5.0, 5.4, 3.4, 1.6))
  d <- density(rr[rr < 3], from = 0, to = 3)
  ymax <- max(d$y)
  if (SHOW_DROP97) {
    dd <- density(rrd[rrd < 3], from = 0, to = 3)
    ymax <- max(ymax, dd$y)
  }
  plot(NA, xlim = c(0, 3), ylim = c(0, ymax * 1.10), xaxt = "n", yaxt = "n",
       bty = "n", xlab = "Modern recruitment / 1990s recruitment", ylab = "Posterior density")
  k <- d$x <= 1
  polygon(c(d$x[k], 1, 0), c(d$y[k], 0, 0), col = fade(BLU, .34), border = NA)
  polygon(c(d$x[!k], 3, 1), c(d$y[!k], 0, 0), col = fade(GRY, .22), border = NA)
  lines(d, col = INK, lwd = 3)
  if (SHOW_DROP97) lines(dd, col = MID, lwd = 2.6, lty = 5)
  segments(1, 0, 1, ymax * 0.70, col = INK, lwd = 2.5, lty = 2)
  text(1, ymax * 0.74, "no change", cex = 1.15, col = INK, font = 2)
  axis(1, 0:3, lwd = 0, lwd.ticks = 1); axis(2, labels = FALSE, lwd = 0)
  if (SHOW_DROP97) {
    segments(1.30, ymax * 1.02, 1.42, ymax * 1.02, col = INK, lwd = 3)
    segments(1.30, ymax * 0.70, 1.42, ymax * 0.70, col = MID, lwd = 2.6, lty = 5)
  }
  text(1.45, ymax * 1.02, sprintf("P(decline) = %.2f", mean(rr < 1)),
       adj = 0, cex = 1.45, font = 2, col = BLU)
  text(1.45, ymax * 0.88,
       sprintf("ratio %.2f (%.2f\u2013%.2f)", mean(rr),
               quantile(rr, .025), quantile(rr, .975)),
       adj = 0, cex = 1.2, col = MID)
  if (SHOW_DROP97) {
    text(1.45, ymax * 0.70, sprintf("excl. 1997:  P = %.2f", mean(rrd < 1)),
         adj = 0, cex = 1.2, font = 2, col = MID)
    text(1.45, ymax * 0.58,
         sprintf("ratio %.2f (%.2f\u2013%.2f)", mean(rrd),
                 quantile(rrd, .025), quantile(rrd, .975)),
         adj = 0, cex = 1.2, col = MID)
  }
  mtext("How much did it change?", side = 3, line = 1.4, adj = 0,
        cex = 1.35, font = 2, col = INK)
  mtext("interval spans 1, so the decline is suggestive, not established",
        side = 3, line = 0.1, adj = 0, cex = 1.05, col = MID)
}

W <- 1620; H <- 720                       # logical px; 2.25:1, suits a 16:9 slide
as_png <- function(f, show, bg, scale = 2) {
  png(f, width = W * scale, height = H * scale, res = 100 * scale,
      bg = bg, type = "cairo-png")
  on.exit(invisible(dev.off())); body(show, bg)
}
as_pdf <- function(f, show) {
  cairo_pdf(f, width = W / 100, height = H / 100, bg = "transparent")
  on.exit(invisible(dev.off())); body(show, "transparent")
}

as_png("manuscript/figures/era_recruitment.png", TRUE, DECK_BG)
as_png("manuscript/figures/era_recruitment_simple.png", FALSE, DECK_BG)
as_png("manuscript/figures/era_recruitment_talk.png", TRUE, "transparent", scale = 3)
as_pdf("manuscript/figures/era_recruitment_talk.pdf", TRUE)
as_png("manuscript/figures/era_recruitment_simple_talk.png", FALSE, "transparent", scale = 3)
as_pdf("manuscript/figures/era_recruitment_simple_talk.pdf", FALSE)
cat(sprintf("wrote 6 files (1997 sensitivity overlay: %s)\n", SHOW_DROP97))