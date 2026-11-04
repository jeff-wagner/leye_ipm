# Lesser Yellowlegs Integrated Population Model

Demographic analysis of a breeding population of Lesser Yellowlegs (*Tringa flavipes*)
in the Anchorage, Alaska area. The core product is an integrated population model (IPM)
fitted in [NIMBLE](https://r-nimble.org), which combines volunteer wetland counts,
mark–resight encounter histories, and nest productivity data into a single likelihood
so that all three sources inform a shared set of demographic rates.

The repository also contains standalone abundance-trend models that predate the IPM and
are retained for comparison.

---

## Quick start

Everything expects the **repository root** as the working directory — scripts use
relative paths such as `data/…`.

```bash
Rscript scripts/LEYE_juv90_module.R
```

```bash
Rscript scripts/LEYE_IPM.R
```

The IPM takes roughly 10–20 minutes to fit serially. For production runs use the
parallel wrapper, which runs one chain per core:

```bash
Rscript scripts/LEYE_IPM_run_parallel.R
```

Either route caches the posterior to `LEYE_IPM_samples.rds`. Then build the
shareable report:

```bash
quarto render reports/leye_ipm_report.qmd
```

**Requirements:** R (developed against 4.6.x) with `nimble`, `coda`, `MCMCvis`,
`tidyverse`, and `readxl`. The data-preparation scripts additionally need
dplyr ≥ 1.2 (they use `filter_out()`). The report additionally needs Quarto plus
`gt`, `ggplot2`, `patchwork` and `here`.

---

## Reports

**`reports/leye_ipm_report.qmd` is the artefact to share with collaborators.** It
reads the cached posterior — it does *not* refit — so rendering takes seconds.

It produces a single self-contained HTML file (`embed-resources: true`), so the
output can be emailed or dropped in a shared folder and will open anywhere with
no supporting files. Figures are embedded, code is folded away by default, and
the interpretation and caveats travel with the numbers.

By default it picks the most recent `LEYE_IPM_samples*.rds` in the project root.
To point it at a specific fit:

```bash
quarto render reports/leye_ipm_report.qmd -P samples_path:LEYE_IPM_samples.rds
```

The rendered HTML **is** tracked, so collaborators who don't run R can open the
current results directly from the repo. Because each render rewrites the whole
self-contained file, re-render when the results actually change rather than
routinely, to keep history from filling with near-duplicate copies.

`summarize_ipm()` in `LEYE_IPM.R` remains as a **console** quick-look for
checking a fit at the terminal. It no longer writes HTML: two HTML generators
meant two things to keep in sync, and the in-script one could not embed its own
figures.

---

## Repository layout

### `scripts/` — analysis

| File | Purpose |
|---|---|
| `LEYE_IPM.R` | **The main model.** Data prep, NIMBLE model, fit, summaries and figures. |
| `LEYE_IPM_run_parallel.R` | Runs the IPM with one MCMC chain per core. Sources `LEYE_IPM.R` with `SKIP_RUN <- TRUE` so only the model objects are defined. |
| `LEYE_IPM_sensitivity.R` | Refits the IPM with **one** thing changed, so companion fits stay comparable to the main one by construction. Recruitment-prior settings: `cohortmean`, `drop97`, `vague`. Count-stream settings: `full`, `balanced9` (the Anchorage bowl is the main model, not a setting). `dryrun` validates every setting's indexing without fitting. Writes to deliberately non-canonical filenames. |
| `compare_count_subsets.R` | Compares the main fit with the `full` and `balanced9` count-stream companions — λ, trajectory, vital rates, convergence — and checks whether `site_re` tracks observed wetland richness or is shrunk toward the mean. |
| `LEYE_recruitment_two_era.R` | **Supplies the IPM's recruitment priors.** Fits 1990s and modern chick cohorts jointly with shared detection, and estimates the change between eras. |
| `LEYE_recruitment_drop_cohort.R` | Leave-one-cohort-out sensitivity on the two-era recruitment fit. Inherits priors, detection, adult survival and both likelihoods from `LEYE_recruitment_two_era.R`; only the dropped cohorts' chicks and the cohort index change. Defaults to 1997, the cohort that returned 7 of 11 chicks and sets the 1990s baseline. |
| `LEYE_juv90_module.R` | The earlier 1990s-only recruitment module. Superseded for prior generation, but retained: it holds the validated `dChick`/`dAdult` likelihoods that the two-era model sources. |
| `create_chick_eh_2018_2026.R` | Cleans the modern chick resights into an analysis-ready encounter history. |
| `LEYE_tLTRE.R` | Transient LTRE. Reads a cached posterior; decomposes mean λ across pathways and Var(λ) across drivers, separating demographic stochasticity. |
| `LEYE_tag_effects_CJS.R` | Standalone CJS robustness check on the tag effects, independent of the IPM. Tests transience and time-since-tagging confounds, and prints the exposure, sex, banding-era and mass-at-banding comparisons the manuscript quotes as confounds. Set `SKIP_TAG_RUN <- TRUE` to reprint every number from the cache without a twenty-minute refit. |
| `fig_era_recruitment.R` | Builds the era-comparison figure (F7/F7b) from the two-era posterior: modeled recruitment per era with individual cohorts, and the posterior of the modern/1990s ratio. Emits deck and transparent talk variants, PNG and PDF, from one draw. |
| `fig_study_area.R` | Builds the study-area figure (F1): upper Cook Inlet with all 32 volunteer survey wetlands, an Anchorage-bowl enlargement with wetland IDs, and an Alaska locator. Also writes `study_area_bowl`, a two-panel version (the 27 fitted wetlands, unlabelled, plus the locator) that the report uses, and the numbered site key. |
| `make_gis_basemap.R` | One-off: clips the GSHHG full-resolution shoreline to the study extent and caches it as `data/gis/cook_inlet_basemap.gpkg`, so F1 redraws offline. Only needed if the extent changes. |
| `create_EH_2026.R` | Builds the adult encounter history and covariate table from raw banding and resight files. |
| `create_juv_eh_1995_2000.R` | Builds the 1990s chick encounter history from the real banding totals. Birds never resighted after banding are reconstructed rows — **see the caveat below.** |
| `LEYE_trend*.R` | Pre-IPM abundance-trend models (negative-binomial GLMMs). Independent of the IPM — and note they fit **all 32** wetlands, where the IPM now fits the 27 in the Anchorage bowl, so their trends are not on the same spatial footprint as λ. |

### `reports/` — collaborator-facing output

| File | Purpose |
|---|---|
| `leye_ipm_report.qmd` | Quarto results report. Reads the cached posterior; does not refit. |
| `report.css` | Styling for the report's headline result cards. |

### `manuscript/`

| File | Purpose |
|---|---|
| `PLAN.md` | Working manuscript plan — target journal, claims, blocking items, section outline, figure list. Tick items as they complete. |

### Cached posteriors (all gitignored)

| File | Contents |
|---|---|
| `LEYE_IPM_samples.rds` | **The canonical fit.** The report and tLTRE read this exact name. |
| `LEYE_IPM_samples_vaguepriors.rds` | Companion fit with flat recruitment priors, for the sensitivity section. |
| `LEYE_IPM_samples_cohortmean.rds` | Companion fit with `pi_rec` centred on the cohort-averaged modern recruitment. |
| `LEYE_IPM_samples_drop1997.rds` | Companion fit with `pi_rec` and `kappa` re-derived from the two-era fit with 1997 excluded. |
| `LEYE_recruitment_two_era_drop1997.rds` | The two-era recruitment posterior with 1997 removed. Supplies the prior above, and the sensitivity overlay in the era figure. |
| `LEYE_IPM_samples_full.rds` | Count-stream companion with all 32 wetlands — the wider-area sensitivity, i.e. the count stream as it was before the main model was restricted to the Anchorage bowl. |
| `LEYE_IPM_samples_balanced9.rds` | Counts restricted to the 9 wetlands surveyed in all 13 years, so panel composition is constant. Not nested in the main fit: two of the nine are outside the bowl. |
| `LEYE_IPM_samples_omega_timevarying.rds` | Archived evidence for the rejected time-varying immigration test. |
| `LEYE_tLTRE_results.rds`, `LEYE_tag_effects_CJS.rds` | Outputs of the corresponding scripts. |

The report and tLTRE require the **exact** filename `LEYE_IPM_samples.rds` and fail
loudly if it is missing. They previously globbed `LEYE_IPM_samples*.rds` and took the
newest, which silently matched the companion fits above — a report can render cleanly
from the wrong posterior that way, so the glob was removed deliberately. Point either
at a different fit explicitly (`-P samples_path:…` for the report, the `samples_path`
argument for the tLTRE) rather than by renaming files.

### `data/`

**Derived — regenerated by the scripts above, not edited by hand:**

| File | Contents |
|---|---|
| `LEYE_95_26_EH.csv` | Adult encounter history, 391 birds × 32 annual occasions (1995–2026). |
| `LEYE_95_26_covs.csv` | Individual covariates, **row-aligned with the EH**. |
| `LEYE_juv_EH_1995_2000.csv` | 1990s chick encounter history, 122 birds × 6 occasions (real banding totals, recovered 2026-09-29). |
| `LEYE_chick_EH_2018_2026.csv` | Modern chick encounter history, 133 birds × 9 occasions, built by `create_chick_eh_2018_2026.R`. |

`LEYE_95_26_EH.csv` and `LEYE_95_26_covs.csv` must stay in **lockstep** — same band
numbers, same order. `LEYE_IPM.R` asserts this and stops if it fails. Regenerate both
together with `create_EH_2026.R`; never edit one alone.

**Raw inputs:**

| File | Contents |
|---|---|
| `New Bands 2025.csv` | 2025 banding records, including an `Age` column (Adult/Chick). |
| `Resight 2026.csv` | 2026 field resightings, keyed by flag code. |
| `Alaska banded LEYE_2018-2025.xlsx` | Flag code → band number lookup for tagged birds. |
| `AnnualSucc_ProdData_95_26.csv` | Nest productivity; `SpecChickNumber` = chicks per nest. Non-missing 2018–2026. |
| `chicks_banded_2018_2026.csv` | Chicks banded per year during the study period. |
| `leye_90s_chick_survival.xlsx` | 1990s chick resightings plus the real 1996 brood data. |
| `Modeled_*.csv` | Distance-sampling data. **Not read by any script in this repo.** |

**Spatial:**

| File | Contents |
|---|---|
| `gis/cook_inlet_basemap.gpkg` | Shoreline, lakes and rivers of upper Cook Inlet, clipped from GSHHG full resolution by `make_gis_basemap.R` (180 KB). Cached so F1 redraws without the 142 MB source archive or an internet connection. Cite Wessel & Smith (1996). |
| `bnb_sites.csv` | The 32 survey wetlands with coordinates, mirrored from `Arin LEYE R code/BnB site metadata.xlsx` on every run of `fig_study_area.R` and checked against it when both are present. |

Count data lives in `Arin LEYE R code/data_2026.csv` (935 wetland × week × year
records; 2014–2026, 32 wetlands, 3 weekly visits), and the wetland coordinates in
`BnB site metadata.xlsx` beside it. That folder holds the original collaborator
analysis and is the source of the abundance series.

**The IPM fits 780 of those 935 records, over 27 of the 32 wetlands.** Since
2026-11-04 the count stream is restricted to the Anchorage bowl so that it covers the
same ground as the mark–resight and nest data; the five dropped wetlands lie 37–67 km
out. The list and the reasoning are at the top of `LEYE_IPM.R`, and the wider-area fit
is kept as a sensitivity. See the caveat below for what the restriction does to λ.

---

## The model

Three data sources share common demographic parameters:

1. **Abundance** — volunteer wetland counts, 2014–2026, as a *relative* index.
2. **Survival** — adult CJS encounter histories, 2010–2026, with individual covariates
   for sex (F/M/Unknown) and tag type (none/geolocator/GPS).
3. **Productivity** — chicks per nest, 2018–2026.

**Survival** is modeled on the logit scale with female and untagged as reference levels:

```
logit(phi[i,t]) = lphi_base[t] + bSex_phi[sex_i] + bTag_phi[tag_i]
logit(p[i])     = lp_base      + bTag_p[tag_i]
```

`phi_ad[t] = ilogit(lphi_base[t])` is the population rate driving the process model
(the population is overwhelmingly untagged). Tag effects are constant in time — the
data are too sparse for a tag × year interaction.

**Recruitment** uses a two-year lag. Chicks are absent from the breeding site until they
recruit, and the 1990s data show recruitment split roughly evenly between ages 1 and 2:

```
R_mean[t] = N_ad[t-1] * (f_rate[t-1]/2) * psi1
          + N_ad[t-2] * (f_rate[t-2]/2) * psi2
```

`psi1` and `psi2` are **composite survive-and-return-locally probabilities, not annual
survival rates.** Birds that survive but settle elsewhere are not counted; immigration
from other populations enters separately through `omega`. Do not compare `psi` against
published juvenile survival estimates — they measure different things.

There is no separate juvenile stage in `N_tot`, because volunteer counts cannot
distinguish non-breeding age-1 birds from breeding adults.

**Productivity** is modeled in two parts — whether a nest succeeded, and how many of a
four-egg clutch fledged given success:

```
succ[m]   ~ Bernoulli(s_nest[year])
chicks[m] ~ Binomial(clutch = 4, q_brood * succ[m])
f_rate[t] <- s_nest[t] * clutch * q_brood
```

This replaced a single Poisson on chicks per nest, which **failed its posterior
predictive check outright** (Bayesian p = 0.000). The observed data are 0, 2, 3 or 4
chicks with *no nest producing exactly one* — a four-egg-clutch shorebird either loses
the nest or fledges a near-full brood, and a Poisson cannot produce that shape at any
mean. Note that `var/mean = 1.18`, so a routine overdispersion check would have passed
it; the misfit was in shape, not variance. `f_rate` keeps its meaning (mean chicks per
nest), so the recruitment equation and the tLTRE were unaffected.

### Year-axis bookkeeping

The three data sources end in different years, which the code handles with pooled
"overflow" parameters rather than by truncating data:

- `phi_pre` — pooled adult survival before 2014 (before the count series starts).
- `phi_post` — pooled adult survival for the 2025→2026 interval, which falls *after*
  the count series ends. Counts currently stop at 2025.
- `f_rate[T+1]` — the 2026 fecundity estimate, which has no corresponding count year and
  so does not feed the process model.

When the count series is extended, these pooled slots shrink accordingly.

### Model checking

**Goodness of fit.** Each data stream carries its own Freeman–Tukey posterior predictive
check (`FT_*_obs` / `FT_*_rep`, monitored by default). Per-stream rather than pooled,
because an IPM can fit one source well and another badly and a single statistic hides
that. Current Bayesian p-values: counts 0.31, productivity 0.50, mark–resight 0.48.

**Prior sensitivity.** The recruitment prior hyperparameters are passed as *constants*
(`a_pi`, `b_pi`, `a_kap`, `b_kap`), so refitting under vague priors is a config change
rather than a forked copy of the model:

```r
cv <- constants; cv$a_pi <- cv$b_pi <- cv$a_kap <- cv$b_kap <- 1
```

Result, and **this changed when the modern recruitment prior replaced the 1990s one**:

| | Informative | Vague |
|---|---|---|
| `pi_rec` | 0.050 (0.009–0.120) | **0.105** (0.004–0.258) |
| `kappa` | 0.481 (0.252–0.718) | 0.491 (0.024–0.975) |

**The recruitment level is now prior-sensitive.** Under the 1990s-derived prior the
informative and vague fits agreed almost exactly (0.107 vs 0.108), so the level was
safely data-driven. Under the modern prior they differ two-fold, because the count data
and the modern chick resights **disagree**: the counts want recruitment near 0.105, the
chick data says 0.052, and the tighter modern prior wins.

That disagreement is a finding, not a defect, and it should be reported rather than
smoothed over. Three readings are possible and the analysis cannot yet separate them:
the counts cannot distinguish recruitment from immigration, so "recruitment" inferred
from counts absorbs any unmodeled inflow; the assumption that modern chick detection
matches the 1990s may be too generous, which would bias `psi` down; or local recruitment
genuinely is low and immigration genuinely is high.

The **timing** remains prior-driven — `kappa`'s interval still widens from (0.25–0.72)
to (0.02–0.97) with its centre unmoved.

---

## Interpreting output

**Abundance is on a relative scale.** `logB0` and the absolute level of `N_tot` are
non-identifiable — only their sum enters the count linear predictor. Expect `logB0` to
show poor convergence (Rhat well above 1.05); this is structural, not a bug, and does
not affect λ. **Read trajectory shape and λ, never absolute abundance.** Its effective
sample size is ~100 out of 21,000 draws at the default chain length, in every fit, so
its Rhat is itself noisy there and has landed anywhere from 1.01 to 1.19 on luck alone;
tripling the chain fixes it. A slice sampler does not — the ridge needs `logB0` and the
whole discrete `N_ad` trajectory to move together.

**There are two averages of λ, and only one is the growth rate.** Use the **geometric**
mean across years, `exp(rowMeans(log(lambda)))` within each draw: it equals
(`N_tot[T]`/`N_tot[1]`)^(1/(T−1)) exactly, and currently gives 0.930 (0.884–0.977),
P(λ<1) = 0.997. The **arithmetic** mean is about 0.023 higher (0.953) because of
year-to-year variation in λ, which compounds to a 14-point difference over twelve years,
so it must never be quoted as a growth rate — it was, until 2026-11-03. It is correct in
exactly one place: `tltre_level()`, where λ = survival + recruitment + immigration holds
exactly each year and only the arithmetic mean preserves that identity, so the pathway
shares sum to 100%. Do not "fix" that one to match.

**Recruitment and immigration are confounded.** Both `psi` and `omega` add birds, and
relative-scale counts cannot fully separate them. This is why recruitment carries an
informative prior from the two-era module: without it, the split is essentially
unidentified. Even with it, the share of annual gains attributable to local recruitment
has a wide credible interval. Report it as such.

Because the two trade off, sharpening one moves the other. Replacing the 1990s-derived
recruitment prior with the modern-era estimate halved recruitment and raised immigration
to compensate — the share of annual gains from local recruitment fell from **56% to 31%**
while λ barely moved (0.966 → 0.965, arithmetic-mean λ as reported at the time; see the
note on the two λ summaries below). The trajectory is robust; the *attribution* is not,
and it depends on which era the recruitment prior comes from.

**The counts carry almost no information about recruitment.** Checked against the fit
rather than assumed, 2026-11-04 — this paragraph previously said `pi_rec` "moves
appreciably from its prior", which is not true of any fit in this repository. On the
current fit `pi_rec` shifts −0.07 prior SDs and the data remove ~10% of its prior
spread; `kappa` (the age-1 versus age-2 split) shifts −0.02 prior SDs and its interval
is unchanged (0.464 → 0.466 wide). Read the recruitment estimates as the chick
analysis's answer, lightly updated — not as something the counts established — and say
so in any write-up. The report quotes these figures inline from the fit so they cannot
go stale again.

---

## Caveats

### The 1990s chick encounter history is partly reconstructed

Cohort **totals are real measurements.** The number of chicks banded in each natal year
was recovered from the original banding records on 2026-09-29 and is now a documented
constant in `create_juv_eh_1995_2000.R`, guarded so the output file cannot disagree
with it:

```r
CHICKS_BANDED <- c(`1995` = 37, `1996` = 34, `1997` = 11, `1998` = 30, `1999` = 10)
```

This replaced an earlier version that drew the 1997 and 1998 denominators from
`rpois(4, 0.85 * num_broods96)` and omitted 1999 altogether. Refitting on the real
numbers left every conclusion unchanged, which is itself a robustness result.

Individual **rows** are still partly reconstructed. Only 16 of the 122 chicks have a
post-banding record in `leye_90s_chick_survival.xlsx`. The other 106 are
banded-but-never-resighted birds, written as a 1 in the natal year and 0 elsewhere.
That they exist is a measurement — the banding totals say how many there were — but
they carry no identity beyond natal year.

What this does and does not limit:

- **Cohort denominators are measurements, not assumptions.** Recruitment scales
  inversely with them, so before 2026-09-29 this was the largest single source of
  uncertainty in the recruitment estimates. It no longer is.
- **The rebuild is deterministic and reproduces the committed data.** Re-running the
  script twice gives byte-identical output, and that output differs from the committed
  CSV only in line endings (the committed file is CRLF, a fresh run writes LF) — every
  parsed value, including row order, is identical. The former warning to treat the CSV
  as the data of record and regenerate only deliberately no longer applies.
- **Brood is unassignable for all but one cohort.** Recorded brood IDs survive only for
  1996 — 12 broods across 34 chicks, in the `broods_1996` sheet. The 106 reconstructed
  rows cannot be grouped at all, so no model in this repo can fit a brood-level random
  effect, and chicks within a cohort are treated as conditionally independent. The
  1990s nicknames (`Firehouse Blue`, `Connors' Kid`) are brood-site labels and do
  identify siblings among the 16 named birds, but not among the rest.
- **Do not use the `Sex` column of the reconstructed rows.** For 1996 the values are
  real sexes taken from the brood records, but *which* reconstructed row receives which
  value is arbitrary; every other cohort is recorded as unknown. No script in this repo
  uses sex as a covariate.
- **Shared detection between 1990s chicks and adults is confirmed, not inferred.** The
  juvenile module fits one detection process across chicks and the 178 contemporaneous
  adults (`scripts/LEYE_recruitment_two_era.R`), which is what turns five modern
  resights out of 133 from a bound into an estimate. That is only legitimate if both
  groups were read the same way. Until 2026-11-04 the justification was circumstantial —
  nicknames like "Firehouse Blue" look like colour-mark combinations — and the
  alternative, metal bands only, would have meant lower true chick detection and an
  underestimated `psi`. **J. Wagner confirmed from the field records on 2026-11-04 that
  the 1990s chicks did carry colour marks.** The assumption stands as written, so
  nothing was refitted; what changed is that it no longer rests on an inference from
  birds' nicknames.

The remaining limitation on the 1990s half of the era comparison is no longer the
denominators but their **concentration**: the 1997 cohort is 11 chicks with 7 returning,
and it sets the baseline the era contrast is measured against. See the recruitment
caveat in the report, and `scripts/LEYE_recruitment_drop_cohort.R` for the
leave-one-cohort-out refit that quantifies it.

### Other open items

- **The three data streams now cover the same ground — this was a deliberate change.**
  The volunteer survey covers 32 wetlands spanning 67 km, but every marked bird and
  every monitored nest comes from the Anchorage bowl. Until 2026-11-04 the IPM fitted
  all 32, so λ described a wider population than the rates meant to explain it. The
  count stream is now restricted to the 27 bowl wetlands (`IPM_COUNT_WETLANDS` in
  `LEYE_IPM.R`), dropping both Palmer Hay Flats accesses, Eklutna, Big Lake and Lynn
  Fuller's Ponds — 37–67 km out, 22.6% of all birds counted on 16.6% of surveys.
  **This was not a way of removing an inconvenient signal:** the restriction makes the
  decline *steeper*, λ 0.939 → 0.930, with P(λ<1) = 0.997 either way. The wider-area
  fit is retained as the `full` sensitivity.
- **The survey panel is unstable, and that was tested too.** It grew from 13 wetlands
  in 2014 to 32 by 2021, only 9 were surveyed throughout, and the ones added late are
  the rich ones (correlation between a site's mean count and its mean survey year,
  +0.50). The site-level random effect handles it: `site_re` tracks observed richness
  with r = 0.98 and slope 0.83 across the 29 wetlands with any detection, so panel
  growth is not leaking abundance into the trend. A nine-wetland balanced panel
  (`balanced9`) flattens λ to 0.950 with the interval crossing 1 — but on a third of
  the data, with a 52% wider interval, so it bounds how much the claim leans on panel
  size rather than contradicting it.
- **Immigration is constant by design, and this was tested.** A time-varying `omega[t]`
  (year random effect) was fitted and rejected: `sd_omega`'s posterior barely moved off
  its prior, annual rates spanned only 1.27×, and in the tLTRE the driver-explained share
  of Var(λ) rose 20.9% → 28.2% while the demographic term fell by the same amount — the
  year effect was absorbing Poisson noise, not detecting signal. It also shifted `omega`
  down and `pi_rec` up without evidence, since the two enter λ additively. The cause is
  arithmetic: ~3 immigrants a year, ~10 usable growth increments. See the comment block
  at the `omega` prior in `LEYE_IPM.R`. Resolving this needs data separating immigrants
  from local recruits, not a richer model.
- **Study-period chick resightings are not used.** Only six unflagged birds were resighted
  across 2018–2026, and detection for unflagged birds is not identifiable from a single
  recapture. Including them would mostly launder a prior assumption into the results.
- **Counts end at 2025** while survival and productivity run to 2026.
- **Older cached samples may be a plain list, not an `mcmc.list`.** The parallel
  wrapper's `lapply()` used to strip the class, so `as.matrix()` on a reloaded cache
  silently returned a 3 × 1 matrix instead of the full draws — without erroring. Fixed at
  the source, and the report and tLTRE both re-wrap defensively, but any cache saved
  before that fix needs `coda::as.mcmc.list()` on load.

---

## Regenerating derived data

```bash
Rscript scripts/create_EH_2026.R
```

Rewrites `LEYE_95_26_EH.csv` and `LEYE_95_26_covs.csv` together. Chicks are filtered out
at source so they never enter the adult survival analysis — a bird banded as a chick is
not an adult at risk, and including it biases adult survival downward.

```bash
Rscript scripts/create_juv_eh_1995_2000.R
```

Rewrites `LEYE_juv_EH_1995_2000.csv`. Deterministic: cohort totals are fixed constants
and the one remaining `slice_sample()` is seeded, so a rerun reproduces the committed
data exactly (bar CRLF versus LF line endings). See the caveat above for what the file
does and does not contain.

---

## Validation

`LEYE_juv90_module.R` uses custom NIMBLE distributions that marginalise over latent
recruitment age. Custom likelihoods fail silently when wrong, so they are checked two
ways: against brute-force enumeration over all latent states, and by confirming they sum
to exactly 1 over the full space of possible encounter histories. The second check found
a real error during development that the first did not — agreement with a reference
implementation means nothing if both share the same mistake. **Re-run both checks after
any change to those functions.**
