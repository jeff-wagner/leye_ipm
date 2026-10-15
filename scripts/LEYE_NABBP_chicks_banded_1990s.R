library(tidyverse)
library(terra)
library(tmap)
library(gt)

nabbp_grp14 <- read_csv("data/NABBP_2025_grp_14.csv")
species_lookup <- read_csv("data/species.csv")
age_lookup <- read_csv("data/age.csv")

leye_id <- species_lookup |>
  filter(ALPHA_CODE == "LEYE") |>
  select(SPECIES_ID) |>
  as.numeric()

leye_ak_juv_90s <- nabbp_grp14 |>
  filter(species_id %in% leye_id & iso_subdivision == "US-AK") |>
  filter(event_year >= 1990 & event_year < 2000) |>
  filter(age_code >= 2 & age_code <= 5)

leye_ak_juv_90s |>
  group_by(event_year) |>
  summarize(n = n())

# Map
sp_leye <- vect(
  leye_ak_juv_90s,
  geom = c("lon_dd", "lat_dd"),
  crs = "EPSG:4326"
)
tmap_mode("view")
tm_shape(sp_leye) +
  tm_dots()

# Anchorage bbox
anc_bbox <- ext(-150, -149.74, 61.07, 61.26)
sp_leye_anc <- crop(sp_leye, anc_bbox)

tm_shape(sp_leye_anc) +
  tm_dots()

sp_leye_anc |>
  as.tibble() |>
  group_by(event_year) |>
  summarize(n = n()) |>
  gt()
