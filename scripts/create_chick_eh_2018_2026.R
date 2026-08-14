# =============================================================================
# Modern-era chick encounter history, 2018-2026
#
# Cleans data/chick_resights_2019_2026_eh.csv into an analysis-ready encounter
# history for the recruitment module. Corrections are applied in code rather
# than by hand so the provenance of the analysed data stays readable.
#
# SOURCE OF TRUTH: this file defines the cohorts. data/chicks_banded_2018_2026.csv
# is retained only as a cross-check and is NOT used to size cohorts -- the two
# disagree for 2018 (23 birds here vs 9 there) and 2025 (5 vs 6, the latter
# being a chick that was never banded). The resight file is authoritative.
# =============================================================================

library(tidyverse)

# band_number is read as CHARACTER, not numeric: two birds are identified by a
# colour-band combination code (dguldgll1, dguldgll2) rather than a metal band
# number, which is consistent with detection by scope read. Treating these IDs
# as numeric would coerce them to NA and silently lose two birds.
raw <- read_csv("data/chick_resights_2019_2026_eh.csv",
  show_col_types = FALSE,
  col_types = cols(band_number = col_character())
) |>
  mutate(band_number = trimws(band_number))
YRS <- 2018:2026
ycols <- as.character(YRS)

# ---- 1. one row per bird ----------------------------------------------------
# An earlier version of this file carried band 129231993 twice, with the two
# rows disagreeing about whether it was resighted. Stop rather than silently
# picking if that ever recurs.
dups <- raw |>
  count(band_number) |>
  filter(n > 1)
if (nrow(dups)) {
  stop(
    "Duplicate band numbers in the resight file: ",
    paste(dups$band_number, collapse = ", "),
    ". Each bird needs a single reconciled encounter history."
  )
}

# ---- 2. band-number typo ----------------------------------------------------
# New Bands 2025 records this chick as 1422220944 (10 digits). Every other band
# in the series is a 9-digit 14222xxxx, so the long form carries a doubled "2".
# The matching adult-side error (1422220945) is fixed in create_EH_2026.R.
BAND_FIX <- c(`1422220944` = "142220944")

# ---- 3. band 129232054: a confirmed local recruit ---------------------------
# Banded as a chick in 2023, then RECAPTURED as a breeding adult in 2025, which
# is why it also appeared in New Bands 2025. The two records disagreed about
# 2025 (chick file 0, adult record 1); the recapture is a real encounter, so the
# reconciled history is the union. This is the only recapture in the data set --
# every other detection was a scope read of colour bands.
#
# create_EH_2026.R drops this bird from the 2025 new-bandings so it is not
# counted twice, and so its first capture is not recorded four years late.
RECRUIT_BAND <- "129232054"
RECRUIT_ADD_YEARS <- c(2025)

eh <- raw |>
  mutate(band_number = coalesce(BAND_FIX[as.character(band_number)], band_number)) |>
  select(band_number, year_banded, all_of(ycols))

for (y in RECRUIT_ADD_YEARS) {
  eh[eh$band_number == RECRUIT_BAND, as.character(y)] <- 1
}

# occasions before a bird was banded are not at risk
for (y in YRS) {
  eh[[as.character(y)]][eh$year_banded > y] <- NA_integer_
}

# ---- 4. cross-check against the banding tally (informational only) ----------
cb <- read_csv("data/chicks_banded_2018_2026.csv", show_col_types = FALSE)
cmp <- eh |>
  count(year_banded, name = "cohort") |>
  full_join(cb, by = c("year_banded" = "Year")) |>
  mutate(difference = cohort - n_chicks) |>
  arrange(year_banded)
cat("=== cohort sizes (this file) vs chicks_banded tally ===\n")
print(as.data.frame(cmp), row.names = FALSE)
cat("Cohorts come from THIS file; differences above are recorded, not resolved.\n")

# ---- 5. summarise -----------------------------------------------------------
M <- as.matrix(eh[, ycols])
bi <- match(eh$year_banded, YRS)
at_risk <- bi < length(YRS) # a bird banded in the final year has no interval
ever <- vapply(seq_len(nrow(eh)), function(i) {
  if (!at_risk[i]) return(FALSE)
  any(M[i, (bi[i] + 1):length(YRS)] == 1, na.rm = TRUE)
}, logical(1))

cat(sprintf(
  "\nchicks %d | at risk %d | ever resighted %d (%.1f%% of at-risk)\n",
  nrow(eh), sum(at_risk), sum(ever), 100 * sum(ever) / sum(at_risk)
))

first_ret <- vapply(which(ever), function(i) {
  which(M[i, (bi[i] + 1):length(YRS)] == 1)[1]
}, numeric(1))
cat("age at first return:",
  paste(names(table(first_ret)), table(first_ret), sep = " -> ", collapse = ",  "),
  "\n"
)
cat("(1990s module, for comparison: 1 -> 6,  2 -> 6,  3 -> 1)\n")

cat("\nresighted individuals:\n")
print(as.data.frame(eh[ever, ]), row.names = FALSE)

write_csv(eh, "data/LEYE_chick_EH_2018_2026.csv")
cat("\nWrote: data/LEYE_chick_EH_2018_2026.csv\n")
