# ========================================================
# Title: LEYE data cleaning
# Author: Jeff Wagner
# Date: 2026-07-31
# Description: Formatting LEYE banding and resight data into encounter history for analysis
# ========================================================

# ---- Set Up ----
library(tidyverse)


# Encounter History ------------------------------------------------------

newBands25 <- read_csv("data/New Bands 2025.csv") |>
  filter(Age != "Chick")
resight26 <- read_csv("data/Resight 2026.csv")
eh_95_25 <- read_csv("data/LEYE_95_25_EH.csv")
banded_18_25 <- readxl::read_excel(
  "data/Alaska banded LEYE_2018-2025.xlsx",
  sheet = "Anchorage tagged birds 2018-25"
)

# Add new bands from 2025
eh_95_26 <- eh_95_25 |>
  mutate(Year = NA) |>
  add_row(
    Band_Number = newBands25$`Band number`,
    Study_Region = newBands25$`Study area`,
    Year = 2025
  ) |>
  filter_out(is.na(Band_Number)) |> # remove chick that wasn't banded
  mutate(
    yr_25 = case_when(Year == 2025 ~ 1, .default = yr_25),
    across(yr_95:yr_24, ~ if_else(Year == 2025, 0L, .x, missing = .x)),
    Study_Region = replace_values(
      Study_Region,
      "Clitheroe" ~ "CLITHEROE",
      "CLITHROE" ~ "CLITHEROE",
      "Oceanview" ~ "OCEAN",
      "Coastal Refuge" ~ "COASTAL"
    )
  ) |>
  select(-Year)

# Add resights
flags26 <- resight26 |>
  filter(!(`Flag code` %in% c("[none]", "[NONE]", "UNK", "-", "KH or XH"))) |> # remove missing or uncertain flag codes
  left_join(
    banded_18_25 |> select(`Flag Code`, Year, `Band Number`),
    by = c("Flag code" = "Flag Code")
  )

eh_95_26 <- eh_95_26 |>
  mutate(
    yr_26 = case_when(
      Band_Number %in% flags26$`Band Number` ~ 1,
      .default = 0
    )
  ) |>
  arrange(Band_Number)

# QC: Check consistency with 2025 EH
mismatch_matrix <- eh_95_25 != eh_95_26[1:381, 1:33]
which(mismatch_matrix, arr.ind = TRUE) # just mistmatches in study area (2nd column) due to recoding, not a concern

write_csv(eh_95_26, "data/LEYE_95_26_EH.csv")

# Covariates -------------------------------------------------------------
covs25 <- read_csv("data/LEYE_95_25_covs.csv")

covs26 <- covs25 |>
  add_row(
    Band_Number = newBands25$`Band number`,
    Study_Region = newBands25$`Study area`,
    Geo = 1,
    GPS = 1,
    Tag_type = 1,
    Carried_tag = 1,
    init_mass = newBands25$`Mass (g)`,
    mean_wing = newBands25$`Wing length (mm)`,
    mean_tot_head = newBands25$`Total head`,
    mean_tarsus = newBands25$`Diagonal tarsus`,
    Otter_18 = 1,
    Resight_effort = 2,
    Sex.f = NA,
    Sex = newBands25$Sex,
    yr_banded = 2025,
    f = 31,
    h = 31
  ) |>
  filter(!is.na(Band_Number)) |> # remove chick that wasn't banded
  mutate(
    Study_Region = replace_values(
      Study_Region,
      "Clitheroe" ~ "CLITHEROE",
      "CLITHROE" ~ "CLITHEROE",
      "Oceanview" ~ "OCEAN",
      "Coastal Refuge" ~ "COASTAL"
    )
  ) |>
  arrange(Band_Number)

write_csv(covs26, "data/LEYE_95_26_covs.csv")
