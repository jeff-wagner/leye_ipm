# =============================================================================
# Study-area basemap -- one-off extraction, cached in the repository
#
# Builds data/gis/cook_inlet_basemap.gpkg: the shoreline, lakes and rivers of
# the Anchorage / upper Cook Inlet area, clipped from GSHHG full resolution.
# scripts/fig_study_area.R reads only that file, so the figure regenerates
# offline and nobody needs the 142 MB source archive to redraw F1.
#
# WHY GSHHG AND NOT NATURAL EARTH
# rnaturalearth ships with this project's dependencies and would need no
# download, but at 1:50 m its coastline has 37 vertices inside this bounding
# box -- Knik Arm and Turnagain Arm, the two features that make the map legible
# as Anchorage, are barely drawn. GSHHG full resolution is ~100 m and renders
# both correctly. The 1:10 m Natural Earth tier sits in between and is still
# too coarse at a 50 km extent.
#
# SOURCE
#   Wessel, P. & Smith, W.H.F. (1996) A global, self-consistent, hierarchical,
#   high-resolution shoreline database. J. Geophys. Res. 101(B4), 8741-8743.
#   gshhg-shp-2.3.7.zip from https://www.soest.hawaii.edu/pwessel/gshhg/
#   Released under the LGPL; cite the paper in the figure caption.
#
# USAGE -- only needed if the cached .gpkg is rebuilt or the extent changes.
#   1. Download and unzip gshhg-shp-2.3.7.zip somewhere OUTSIDE this repo.
#   2. Rscript scripts/make_gis_basemap.R <path-to-unzipped-dir>
#      (the directory that contains GSHHS_shp/ and WDBII_shp/)
#
# The archive is third-party data, so it stays out of the working tree and is
# passed in by path rather than sourced from a fixed location.
# =============================================================================

suppressMessages(library(sf))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  stop("Usage: Rscript scripts/make_gis_basemap.R <dir containing GSHHS_shp/>")
}
SRC <- normalizePath(args[1], mustWork = TRUE)
OUT <- "data/gis/cook_inlet_basemap.gpkg"

# Extent. The 32 survey wetlands span 61.06-61.64 N, 150.06-149.02 W; this adds
# enough margin to carry Knik Arm, the head of Turnagain Arm, the Palmer Hay
# Flats and the Susitna side of the inlet, which is the context a reader needs
# to place the sites. Widening it later means rerunning this script, not
# editing the figure.
BB <- c(xmin = -150.55, ymin = 60.85, xmax = -148.70, ymax = 61.85)
box <- st_as_sfc(st_bbox(BB, crs = 4326))

# wkt_filter pushes the clip into GDAL, so the 161 MB L1 shapefile is never
# read into memory in full.
grab <- function(rel, kind) {
  p <- file.path(SRC, rel)
  if (!file.exists(p)) {
    message("  skipped (not extracted): ", rel)
    return(NULL)
  }
  g <- st_read(p, wkt_filter = st_as_text(box), quiet = TRUE)
  if (!nrow(g)) {
    message("  nothing inside the extent: ", rel)
    return(NULL)
  }
  g <- st_make_valid(g)
  g <- suppressWarnings(st_intersection(g, box))
  g <- g[!st_is_empty(st_geometry(g)), ]
  if (!nrow(g)) return(NULL)
  out <- st_sf(kind = kind, geometry = st_geometry(g))
  message(sprintf("  %-26s %3d features", rel, nrow(out)))
  out
}

message("Clipping GSHHG to ", paste(BB, collapse = " "), " ...")
layers <- list(
  land = grab("GSHHS_shp/f/GSHHS_f_L1.shp", "land"),
  lake = grab("GSHHS_shp/f/GSHHS_f_L2.shp", "lake"),
  island_in_lake = grab("GSHHS_shp/f/GSHHS_f_L3.shp", "island_in_lake"),
  river = do.call(rbind, Filter(Negate(is.null), lapply(
    sprintf("WDBII_shp/f/WDBII_river_f_L%02d.shp", 1:5), grab, kind = "river"
  )))
)
layers <- Filter(Negate(is.null), layers)
stopifnot("land" %in% names(layers))

if (file.exists(OUT)) unlink(OUT)
for (nm in names(layers)) {
  st_write(layers[[nm]], OUT, layer = nm, quiet = TRUE,
           append = file.exists(OUT))
}
message(sprintf(
  "\nWrote %s (%.0f KB): layers %s",
  OUT, file.size(OUT) / 1024, paste(names(layers), collapse = ", ")
))
