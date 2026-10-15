# ========================================================
# Title: Build 1990s LEYE banded chick encounter history
# Author: Jeff Wagner
# Date: 2026-08-04 (rewritten 2026-09-29)
# Description: Assembles the 1995-1999 chick cohorts from the real banding
#   records: the survival sheet holds every chick with a post-banding history,
#   and each cohort is completed to its true size with never-resighted birds.
# ========================================================

# ---- Set Up ----
library(tidyverse)

# Real chicks banded per natal year, recovered from the original banding
# records on 2026-09-29. These replace the simulated denominators this script
# used to draw with rpois(4, 0.85 * num_broods96); that branch is gone, and so
# is the warning that re-running produced a different, equally arbitrary
# cohort structure. The numbers are now fixed, so the file is reproducible.
#
# 1996 is corroborated independently by the broods_1996 sheet (34 chicks).
# 1999 is a cohort this file previously omitted altogether.
#
# NOTE ON 1997: 11 chicks banded, 8 of them carry a post-banding history and 7
# were resighted. A 64% local return rate is six to sixteen times every other
# cohort in either era, and this one cohort therefore dominates the 1990s
# recruitment signal. 1997 was also a low-effort season for adult marking
# (21 adults, against 52/48/36/21 in the surrounding years), so the small
# denominator is real. Treat psi from the 1990s as resting heavily on it.
CHICKS_BANDED <- c(
  `1995` = 37,
  `1996` = 34,
  `1997` = 11,
  `1998` = 30,
  `1999` = 10
)

# The only randomness left is which real 1996 sexes go to which padded 1996
# bird; cohort sizes no longer depend on the seed.
set.seed(20260804)

eh <- readxl::read_xlsx(
  "data/leye_90s_chick_survival.xlsx",
  sheet = "leye_90s_chick_survival"
)
broods96 <- readxl::read_xlsx(
  "data/leye_90s_chick_survival.xlsx",
  sheet = "broods_1996"
)

# `1998`:`2000` come in as character (mixed "0"/"1"/"NA"); coerce so new rows
# bind cleanly. NA means "not monitored", and is kept -- the downstream module
# drops a bird whose banding occasion is its last monitored one.
eh_clean <- eh |>
  mutate(across(`1998`:`2000`, as.numeric))

band_years <- as.integer(names(CHICKS_BANDED))
occasions <- c(1995, 1996, 1997, 1998, 1999, 2000)

# How many rows each cohort already has, and how many must be added.
have <- eh_clean |>
  count(`Year banded`, name = "n_records") |>
  complete(`Year banded` = band_years, fill = list(n_records = 0L))
need <- CHICKS_BANDED[as.character(have$`Year banded`)] - have$n_records
names(need) <- have$`Year banded`

if (any(need < 0)) {
  stop(
    "More chick records than chicks banded in: ",
    paste(names(need)[need < 0], collapse = ", "),
    ". Either the banding totals or the survival sheet is wrong."
  )
}

# One row per banded-but-never-resighted chick: a 1 in the natal year, 0
# everywhere else, because these birds were never seen again.
make_unresighted_rows <- function(year, sex_vec, name_prefix) {
  n <- length(sex_vec)
  if (n == 0) {
    return(NULL)
  }
  out <- tibble(
    Name = paste0(name_prefix, "_", seq_len(n)),
    `Year banded` = year,
    Loc = "A",
    Sex = sex_vec
  )
  for (y in occasions) out[[as.character(y)]] <- as.numeric(year == y)
  out
}

# 1996 padding takes real sexes from the brood records; other cohorts have no
# sex information for birds never seen again, so they are recorded as unknown
# rather than drawn at random.
sex_1996 <- broods96 |>
  slice_sample(n = need[["1996"]]) |>
  pull(DNA)

new_rows <- bind_rows(
  make_unresighted_rows(1995, rep("U", need[["1995"]]), "Unresighted_1995"),
  make_unresighted_rows(1996, sex_1996, "Unresighted_1996"),
  make_unresighted_rows(1997, rep("U", need[["1997"]]), "Unresighted_1997"),
  make_unresighted_rows(1998, rep("U", need[["1998"]]), "Unresighted_1998"),
  make_unresighted_rows(1999, rep("U", need[["1999"]]), "Unresighted_1999")
)

eh_complete <- bind_rows(eh_clean, new_rows)

# Guard: the file must reproduce the banding totals exactly, or every
# recruitment estimate downstream is scaled wrong.
got <- table(eh_complete$`Year banded`)[as.character(band_years)]
stopifnot(identical(as.integer(got), as.integer(CHICKS_BANDED)))

write_csv(eh_complete, "data/LEYE_juv_EH_1995_2000.csv")

cat(sprintf(
  "Wrote %d chicks across %d cohorts: %s\n",
  nrow(eh_complete), length(band_years),
  paste(names(CHICKS_BANDED), CHICKS_BANDED, sep = "=", collapse = ", ")
))
