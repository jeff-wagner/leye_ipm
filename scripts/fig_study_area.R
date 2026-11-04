# =============================================================================
# F1 -- study area
#
# Run from the PROJECT ROOT:
#     Rscript scripts/fig_study_area.R
#
# Writes manuscript/figures/study_area.png / .pdf and a transparent talk
# variant, all from one body() so the manuscript and the slide cannot drift.
#
# THREE PANELS
#   A  regional -- upper Cook Inlet, all 32 volunteer survey wetlands, with the
#      panel B extent boxed. Bubble area is mean Lesser Yellowlegs per survey.
#   B  Anchorage bowl -- 24 of the 32 wetlands sit inside about 25 km here, so
#      the regional panel cannot carry site numbers without collapsing into
#      noise.
#   C  Alaska locator.
#
# WHAT THE BUBBLES ARE, AND WHAT THEY ARE NOT
# Mean birds per survey is the RAW count, not an abundance estimate: it is not
# corrected for detection, survey effort or wetland area, all of which the trend
# model handles and none of which this figure does. It is here to show that the
# sampling frame spans occupied and unoccupied wetlands -- three sites have
# never produced a detection in thirteen years and are drawn hollow, which is
# information a reader loses if F1 shows undifferentiated dots.
#
# DATA
#   Arin LEYE R code/BnB site metadata.xlsx   wetland names and coordinates
#   Arin LEYE R code/data_2026.csv            the count series the IPM uses
#   data/gis/cook_inlet_basemap.gpkg          shoreline (scripts/make_gis_basemap.R)
#
# The site coordinates are mirrored into data/bnb_sites.csv on every run and
# checked against the spreadsheet when both are present, so the figure has a
# tracked record of the coordinates it drew and the two cannot silently diverge.
#
# NOT ON THIS MAP: the mark-resight and nest-monitoring sites. The encounter
# history carries a `Study_Region` label per bird, but those are field codes
# (OLNO, AMCO, WEPO ...) with no coordinate table in the repository, and
# several of them do not correspond to a surveyed wetland at all. Guessing at
# them from the codes would put invented points on a published figure. Supply a
# site-code -> coordinate table and this script can draw that layer.
# =============================================================================

suppressMessages({
  library(sf)
  library(rnaturalearth)
})

XLSX <- "Arin LEYE R code/BnB site metadata.xlsx"
SITES_CSV <- "data/bnb_sites.csv"
COUNTS <- "Arin LEYE R code/data_2026.csv"
GPKG <- "data/gis/cook_inlet_basemap.gpkg"
OUTDIR <- "manuscript/figures"

## ---------------------------------------------------------------------------
## Data
## ---------------------------------------------------------------------------
read_xlsx_sites <- function() {
  if (!requireNamespace("readxl", quietly = TRUE) || !file.exists(XLSX)) return(NULL)
  d <- suppressMessages(readxl::read_excel(XLSX, sheet = 1))
  d <- data.frame(
    id = as.integer(d[["Wetland ID"]]), name = as.character(d[["Wetland name"]]),
    lat = as.numeric(d[["Latitude"]]), lon = as.numeric(d[["Longitude"]])
  )
  d[!is.na(d$id), ]
}

sites <- read_xlsx_sites()
if (is.null(sites)) {
  if (!file.exists(SITES_CSV)) {
    stop("Need either ", XLSX, " or the mirrored ", SITES_CSV)
  }
  message("Spreadsheet not found; drawing from ", SITES_CSV)
  sites <- read.csv(SITES_CSV)
} else if (file.exists(SITES_CSV)) {
  old <- read.csv(SITES_CSV)
  same <- isTRUE(all.equal(
    old[order(old$id), c("id", "lat", "lon")],
    sites[order(sites$id), c("id", "lat", "lon")],
    check.attributes = FALSE
  ))
  if (!same) message("NOTE: ", SITES_CSV, " disagreed with the spreadsheet and was rewritten.")
}
sites <- sites[order(sites$id), ]
stopifnot(nrow(sites) == 32, !anyNA(sites$lat), !anyNA(sites$lon))

cnt <- read.csv(COUNTS, fileEncoding = "UTF-8-BOM")
agg <- aggregate(count ~ wetland, cnt, function(x) c(mean(x), length(x), sum(x > 0)))
sites$mean_count <- agg$count[match(sites$id, agg$wetland), 1]
sites$n_survey <- agg$count[match(sites$id, agg$wetland), 2]
sites$n_pos <- agg$count[match(sites$id, agg$wetland), 3]
stopifnot(!anyNA(sites$mean_count)) # every wetland must appear in the counts
write.csv(sites[, c("id", "name", "lat", "lon")], SITES_CSV, row.names = FALSE)

YRS <- range(cnt$year)
ZERO <- sites$n_pos == 0

## ---------------------------------------------------------------------------
## Geometry -- everything is drawn in Alaska Albers (EPSG:3338)
## ---------------------------------------------------------------------------
PRJ <- 3338
# st_union on the land: clipping GSHHG to a box splits the mainland into
# fifteen tiled polygons, and drawing their borders draws the seams between
# them as coastline that is not there. Union first, then only real shoreline
# gets a stroke.
base <- lapply(c("land", "lake"), function(l) {
  g <- st_transform(st_read(GPKG, l, quiet = TRUE), PRJ)
  if (l == "land") st_union(g) else st_geometry(g)
})
names(base) <- c("land", "lake")

pt <- st_transform(st_as_sf(sites, coords = c("lon", "lat"), crs = 4326), PRJ)
XY <- st_coordinates(pt)

ll2xy <- function(lon, lat) {
  st_coordinates(st_transform(
    st_sfc(st_multipoint(cbind(lon, lat)), crs = 4326), PRJ
  ))[, 1:2, drop = FALSE]
}
# A lon/lat window is a quadrilateral in Albers, not a rectangle, so take the
# projected bounding box of its densified outline rather than of its corners.
win <- function(lon, lat) {
  e <- seq(0, 1, length.out = 60)
  b <- rbind(
    cbind(lon[1] + diff(lon) * e, lat[1]), cbind(lon[2], lat[1] + diff(lat) * e),
    cbind(lon[1] + diff(lon) * e, lat[2]), cbind(lon[1], lat[1] + diff(lat) * e)
  )
  p <- ll2xy(b[, 1], b[, 2])
  list(x = range(p[, 1]), y = range(p[, 2]))
}
W_REG <- win(c(-150.40, -148.88), c(60.95, 61.78)) # panel A
W_BOWL <- win(c(-150.13, -149.67), c(61.03, 61.33)) # panel B

## ---------------------------------------------------------------------------
## Look -- the palette of the other standalone figures
## ---------------------------------------------------------------------------
VER <- "#D55E00"
BLU <- "#0072B2"
GRY <- "#8A96A0"
DECK_BG <- "#FAF8F4"
INK <- "#1C2B33"
MID <- "#5A6872"
LAND <- "#EDEAE3"
LAND_T <- "#E4E0D6"
WATER <- "#DCE8F1"
SHORE <- "#A9B4BC"
fade <- function(h, a) {
  v <- col2rgb(h)
  rgb(v[1], v[2], v[3], alpha = a * 255, maxColorValue = 255)
}

# Bubble AREA proportional to the mean count, so the eye compares areas and not
# radii; a common scale across panels A and B keeps the two readable together.
RAD <- function(m) 0.9 + 2.6 * sqrt(m / max(sites$mean_count))

# asp = 1 keeps the projection honest, which means the panel almost never has
# the window's exact shape: whichever axis is short gets widened to fill the
# cell. So paint the sea over par("usr"), the window actually drawn, not over
# `w` -- otherwise the figure background shows through as a letterbox and the
# panel looks cropped rather than mapped.
draw_base <- function(w, shore_lwd = 1.1, water = WATER, land = LAND) {
  plot(NA,
    xlim = w$x, ylim = w$y, asp = 1, xaxs = "i", yaxs = "i",
    xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n"
  )
  u <- par("usr")
  rect(u[1], u[3], u[2], u[4], col = water, border = NA)
  plot(base$land, col = land, border = SHORE, lwd = shore_lwd, add = TRUE)
  plot(base$lake, col = water, border = SHORE, lwd = shore_lwd * .8, add = TRUE)
  rect(u[1], u[3], u[2], u[4], col = NA, border = SHORE, lwd = 1.3, xpd = NA)
  invisible(list(x = u[1:2], y = u[3:4]))
}

draw_sites <- function(cx = 1, numbers = FALSE, nudge = NULL) {
  r <- RAD(sites$mean_count) * cx
  points(XY[ZERO, 1], XY[ZERO, 2],
    pch = 21, bg = "white",
    col = GRY, lwd = 1.8, cex = 1.1 * cx
  )
  points(XY[!ZERO, 1], XY[!ZERO, 2],
    pch = 21, bg = fade(BLU, .55),
    col = "white", lwd = 1.3, cex = r[!ZERO]
  )
  if (numbers) {
    off <- if (is.null(nudge)) matrix(0, nrow(XY), 2) else nudge
    text(XY[, 1] + off[, 1], XY[, 2] + off[, 2], sites$id,
      cex = .78 * cx, font = 2, col = INK
    )
  }
}

# True north, not grid north: EPSG:3338 is centred on 154 W, so the meridians
# lean by about 4 degrees here and an arrow drawn straight up would be wrong.
north_arrow <- function(x, y, len) {
  p <- ll2xy(c(-149.6, -149.6), c(61.2, 61.4))
  b <- atan2(p[2, 1] - p[1, 1], p[2, 2] - p[1, 2])
  dx <- sin(b) * len
  dy <- cos(b) * len
  arrows(x, y, x + dx, y + dy, length = .11, lwd = 2.6, col = INK)
  text(x + dx * 1.26, y + dy * 1.26, "N", cex = 1.0, font = 2, col = INK)
}

scale_bar <- function(x, y, km, h) {
  m <- km * 1000
  rect(x, y, x + m / 2, y + h, col = INK, border = INK)
  rect(x + m / 2, y, x + m, y + h, col = "white", border = INK)
  text(c(x, x + m), y - h * 1.1, c("0", paste(km, "km")),
    adj = c(.5, 1), cex = .85, col = INK
  )
}

## ---------------------------------------------------------------------------
## Alaska locator
## ---------------------------------------------------------------------------
ak <- tryCatch(
  {
    us <- ne_countries(
      scale = "medium", country = "united states of america",
      returnclass = "sf"
    )
    st_transform(suppressWarnings(st_crop(
      st_geometry(us),
      st_bbox(c(xmin = -179.9, ymin = 51, xmax = -129, ymax = 72), crs = 4326)
    )), PRJ)
  },
  error = function(e) NULL
)

draw_locator <- function(bg) {
  if (is.null(ak)) {
    plot.new()
    return(invisible())
  }
  b <- st_bbox(ak)
  plot(NA,
    xlim = c(b$xmin, b$xmax), ylim = c(b$ymin, b$ymax), asp = 1,
    xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n"
  )
  plot(ak, col = LAND_T, border = SHORE, lwd = .7, add = TRUE)
  rect(W_REG$x[1], W_REG$y[1], W_REG$x[2], W_REG$y[2],
    border = VER, lwd = 2.4, col = fade(VER, .18)
  )
  text(mean(W_REG$x), W_REG$y[1] - (b$ymax - b$ymin) * .085,
    "study area",
    col = VER, font = 2, cex = .95, adj = c(.5, 1)
  )
  text(b$xmin + (b$xmax - b$xmin) * .16, b$ymin + (b$ymax - b$ymin) * .80,
    "ALASKA",
    col = MID, cex = 1.0, font = 2
  )
}

# Hand-placed label offsets for panel B, in kilometres. Twenty-four wetlands
# inside 25 km will not label themselves; anything automatic here either
# overlaps or wanders far enough from its point to be ambiguous.
NUDGE <- local({
  n <- matrix(0, nrow(XY), 2)
  n[, 2] <- 1.15 # default: just above the point
  set <- function(id, dx, dy) n[match(id, sites$id), ] <<- c(dx, dy)
  set(17, -0.2, 1.35)
  set(13, -1.6, 0.1)
  set(11, 0.4, -1.30)
  set(30, -0.3, 1.30)
  set(21, 1.7, 0.85)
  set(29, -2.0, 0.55)
  set(32, 1.6, 0.55)
  set(7, -1.5, 1.00)
  set(14, -1.6, 0.55)
  set(9, 0.5, -1.35)
  set(12, -0.4, 1.30)
  set(31, 1.5, -0.45)
  set(6, -0.3, 1.30)
  set(1, 1.5, 0.15)
  set(4, -1.4, 0.45)
  set(8, 0.3, 1.30)
  set(16, -1.7, 0.05)
  set(26, -0.4, -1.35)
  set(27, -1.4, -0.45)
  set(28, 1.5, 0.15)
  set(5, -1.4, 0.45)
  set(2, -1.6, 0.15)
  set(19, -1.7, -0.85)
  set(15, 0.4, -1.45)
  set(20, 1.7, 0.25)
  set(25, 1.6, 0.35)
  set(24, 1.5, 0.15)
  n * 1000
})

## ---------------------------------------------------------------------------
## Figure
## ---------------------------------------------------------------------------
body <- function(bg) {
  op <- par(bg = bg, family = "sans", col.axis = MID, col.lab = INK, fg = MID)
  on.exit(par(op))
  land <- if (identical(bg, "transparent")) LAND_T else LAND
  # Column 1 is panel A over both rows; column 2 stacks B above C.
  #
  # DRAWN IN THE ORDER B, C, A -- note the 3 in the layout matrix. Keeping
  # asp = 1 means a panel shows whatever its cell shape makes it show, which is
  # never exactly the window asked for, so the enlarged area in B is wider than
  # W_BOWL. Drawing B first lets panel A outline the ground B ACTUALLY covers
  # instead of the ground it was asked for; otherwise the orange rectangle and
  # the orange panel border would promise the same extent and not deliver it.
  layout(matrix(c(3, 3, 1, 2), 2), widths = c(1.45, 1), heights = c(1.55, 1))

  ## ---- B: Anchorage bowl -------------------------------------------------
  par(mar = c(0.4, 0.4, 2.3, 0.4))
  wB <- draw_base(W_BOWL, shore_lwd = 1.3, land = land)
  draw_sites(cx = 1.15, numbers = TRUE, nudge = NUDGE)
  scale_bar(wB$x[1] + 1800, wB$y[1] + 2000, 5, diff(wB$y) * .018)
  mtext("B   Anchorage bowl",
    side = 3, line = 0.9, adj = 0,
    cex = 1.1, font = 2, col = INK
  )
  mtext("numbers are wetland IDs; names in the site table",
    side = 3, line = -0.1,
    adj = 0, cex = .8, col = MID
  )
  box(col = VER, lwd = 2.2)

  ## ---- C: locator --------------------------------------------------------
  par(mar = c(0.4, 0.4, 2.3, 0.4))
  draw_locator(bg)
  mtext("C   Location",
    side = 3, line = 0.9, adj = 0,
    cex = 1.1, font = 2, col = INK
  )

  ## ---- A: regional -------------------------------------------------------
  par(mar = c(0.4, 0.4, 2.6, 0.4))
  wA <- draw_base(W_REG, land = land)
  rect(wB$x[1], wB$y[1], wB$x[2], wB$y[2], border = VER, lwd = 2.2, lty = 1)
  text(wB$x[2], wB$y[2], "B",
    col = VER, font = 2, cex = 1.15,
    adj = c(-0.35, 0.4)
  )
  draw_sites()

  place <- function(lon, lat, lab, col = MID, cex = .92, font = 3, srt = 0) {
    p <- ll2xy(lon, lat)
    text(p[1], p[2], lab, col = col, cex = cex, font = font, srt = srt)
  }
  # Placed off the panel's own window rather than a lon/lat: the only open
  # water left of the panel B box moves with the aspect fit, and a fixed
  # coordinate here lands under the box or off the edge.
  text(wA$x[1] + diff(wA$x) * .105, wA$y[1] + diff(wA$y) * .115,
    "COOK  INLET",
    col = MID, cex = 1.0, font = 3
  )
  place(-149.90, 61.42, "Knik Arm", srt = 42)
  place(-149.55, 61.03, "Turnagain Arm")
  place(-149.42, 61.55, "Palmer Hay Flats", cex = .88)
  place(-150.14, 61.15, "Fire I.", cex = .78)
  place(-149.82, 61.21, "ANCHORAGE", col = INK, font = 2, cex = .95)
  place(-149.15, 61.15, "CHUGACH  MTS", cex = .88)

  scale_bar(wA$x[1] + 4500, wA$y[1] + 5000, 20, diff(wA$y) * .013)
  north_arrow(wA$x[2] - 9000, wA$y[1] + 6000, diff(wA$y) * .085)
  mtext("A   Upper Cook Inlet",
    side = 3, line = 1.1, adj = 0,
    cex = 1.15, font = 2, col = INK
  )
  mtext(sprintf(
    "32 volunteer-surveyed wetlands, %d–%d   |   bubble area = mean Lesser Yellowlegs per survey   |   hollow = never detected",
    YRS[1], YRS[2]
  ), side = 3, line = 0.05, adj = 0, cex = .82, col = MID)

  # Bubble key, on the empty Susitna-side land at the top left: the only part
  # of panel A with neither a site nor a place name in it.
  key <- c(1, 4, 8)
  kx <- wA$x[1] + diff(wA$x) * .045
  ky <- wA$y[2] - diff(wA$y) * c(.11, .17, .245)
  points(rep(kx, 3), ky,
    pch = 21, bg = fade(BLU, .55), col = "white",
    lwd = 1.3, cex = RAD(key)
  )
  text(rep(kx, 3) + diff(wA$x) * .050, ky, key, cex = .85, col = MID, adj = 0)
  text(kx - diff(wA$x) * .030, ky[1] + diff(wA$y) * .055,
    "birds / survey",
    cex = .85, col = MID, font = 2, adj = 0
  )
}

W <- 1420
H <- 980
as_png <- function(f, bg, scale = 1) {
  png(f,
    width = W * scale, height = H * scale, res = 100 * scale,
    type = "cairo", bg = bg
  )
  on.exit(invisible(dev.off()))
  body(bg)
}
as_pdf <- function(f) {
  cairo_pdf(f, width = W / 100, height = H / 100, bg = "transparent")
  on.exit(invisible(dev.off()))
  body("transparent")
}

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
as_png(file.path(OUTDIR, "study_area.png"), DECK_BG)
as_png(file.path(OUTDIR, "study_area_talk.png"), "transparent", scale = 3)
as_pdf(file.path(OUTDIR, "study_area.pdf"))

# The site key the panel B numbers refer to.
key <- sites[, c("id", "name", "lat", "lon", "n_survey", "n_pos", "mean_count")]
key$mean_count <- round(key$mean_count, 2)
write.csv(key, file.path(OUTDIR, "study_area_sites.csv"), row.names = FALSE)

cat(sprintf(
  "Wrote %s/study_area.{png,pdf}, study_area_talk.png and study_area_sites.csv\n  %d wetlands, %d with no detection in %d-%d\n",
  OUTDIR, nrow(sites), sum(ZERO), YRS[1], YRS[2]
))
