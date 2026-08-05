# ========================================================
# Title: Simulate 1990s LEYE banded chick cohorts
# Author: Jeff Wagner
# Date: 2026-08-04
# Description: This script simulates missing chick banding data from the 1990s period by drawing on available cohort information from 1996.
# ========================================================

# ---- Set Up ----
library(tidyverse)

# Cohort sizes for 1995/1997/1998 are drawn at random (see "Simulate missing
# brood data" below), as are the sexes of never-resighted birds outside 1996.
# Without a seed, every run produces different cohort denominators and so
# different downstream recruitment estimates (psi1/psi2 in LEYE_juv90_module.R
# scale inversely with them). Seed set here so re-runs are reproducible.
#
# NOTE: the committed data/LEYE_juv_EH_1995_2000.csv was generated BEFORE this
# seed was added, so re-running will NOT reproduce that file -- it will produce
# a different, equally valid draw. Re-running therefore invalidates any
# recruitment estimate derived from the committed version.
set.seed(20260804)

eh <- readxl::read_xlsx(
  "data/leye_90s_chick_survival.xlsx",
  sheet = "leye_90s_chick_survival"
)
broods96 <- readxl::read_xlsx(
  "data/leye_90s_chick_survival.xlsx",
  sheet = "broods_1996"
)


# Summary ----------------------------------------------------------------
band_years <- unique(eh$`Year banded`)
resight_years <- as.numeric(colnames(eh[, 5:10]))
num_broods96 <- length(unique(broods96$`Brood name`))
chicks_per_brood96 <- broods96 |>
  group_by(`Brood name`) |>
  summarise(n_chicks = n())


# Simulate missing brood data --------------------------------------------
broods_per_year <- tibble(Year = band_years) |>
  mutate(
    broods = case_when(
      Year == 1996 ~ num_broods96,
      .default = rpois(4, 0.85 * num_broods96) # random poisson draw with arbitrary mean of 0.85*broods banded in 1996
    )
  )

# Simulate chicks per brood using custom categorical distribution based on proportion of chicks/brood in 1996
# Define exact probabilities for 1, 2, 3, and 4 (must sum to 1)
probabilities <- chicks_per_brood96 |>
  group_by(n_chicks) |>
  summarize(prop = n() / num_broods96)
sum(probabilities$prop)

chicks_per_brood <- broods_per_year |>
  rowwise() |>
  mutate(
    chicks = list(
      if (Year == 1996) {
        chicks_per_brood96$n_chicks
      } else {
        sample(1:4, size = broods, replace = TRUE, prob = probabilities$prop)
      }
    )
  ) |>
  ungroup()

# Plot
chicks_per_brood |>
  select(Year, chicks) |>
  unnest_longer(chicks, values_to = "n_chicks") |>
  mutate(brood_id = row_number(), .by = Year) |>
  ggplot(aes(x = n_chicks, y = factor(Year))) +
  geom_jitter(aes(color = Year == 1996), height = 0.15, width = 0.1, size = 2) +
  labs(x = "Chicks per brood", y = "Year", color = "Actual data (1996)") +
  scale_color_manual(values = c("TRUE" = "#4daf3a", "FALSE" = "#3a71e9"))


# Make encounter history for all years --------------------------------------
tot_chicks_per_year <- chicks_per_brood |>
  group_by(Year) |>
  summarize(n_chicks = sum(unlist(chicks)))

resighted_chicks_per_year <- eh |>
  group_by(`Year banded`) |>
  summarize(n_chicks = n())

new_chicks_needed <- tot_chicks_per_year$n_chicks -
  resighted_chicks_per_year$n_chicks

# Add rows for chicks that were banded but never resighted -------------------
# `1998`:`2000` were read in as character (mixed "0"/"1"/"NA" strings);
# coerce to numeric so new rows can be bound cleanly
eh_clean <- eh |>
  mutate(across(`1998`:`2000`, as.numeric))

# One row per never-resighted chick: 0 in all years except a 1 in the
# band year, since these birds were never seen again
make_unresighted_rows <- function(year, sex_vec, name_prefix) {
  n <- length(sex_vec)
  tibble(
    Name = paste0(name_prefix, "_", seq_len(n)),
    `Year banded` = year,
    Loc = "A",
    Sex = sex_vec,
    `1995` = as.numeric(year == 1995),
    `1996` = as.numeric(year == 1996),
    `1997` = as.numeric(year == 1997),
    `1998` = as.numeric(year == 1998),
    `1999` = 0,
    `2000` = 0
  )
}

# Random per-bird sex draws, roughly 50/50 on average
random_sex <- function(n) sample(c("F", "M"), size = n, replace = TRUE)

# 1996 chicks: draw sex from the real broods96 cohort (33 of the 34
# recorded chicks; 1 chick from this cohort is already in eh)
sex_1996 <- broods96 |>
  slice_sample(n = new_chicks_needed[band_years == 1996]) |>
  pull(DNA)

new_rows <- bind_rows(
  make_unresighted_rows(
    1995,
    random_sex(new_chicks_needed[band_years == 1995]),
    "Unresighted_1995"
  ),
  make_unresighted_rows(1996, sex_1996, "Unresighted_1996"),
  make_unresighted_rows(
    1997,
    random_sex(new_chicks_needed[band_years == 1997]),
    "Unresighted_1997"
  ),
  make_unresighted_rows(
    1998,
    random_sex(new_chicks_needed[band_years == 1998]),
    "Unresighted_1998"
  )
)

eh_complete <- bind_rows(eh_clean, new_rows)
write_csv(eh_complete, "data/LEYE_juv_EH_1995_2000.csv")
