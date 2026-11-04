---
title: "Manuscript plan"
subtitle: "Integrated population model, Lesser Yellowlegs, Anchorage"
# gfm_auto_identifiers makes pandoc build heading anchors GitHub's way, so a
# section link copied from GitHub also works in the rendered HTML and vice
# versa. Without it pandoc strips leading digits from headings
# ("1. The argument" -> #the-argument) and the numbered sections get two
# different anchors depending on where the file is being read.
from: markdown+gfm_auto_identifiers
format:
  html:
    toc: true
    toc-depth: 3
    toc-location: left
    toc-title: "Sections"
    toc-expand: true
    embed-resources: true
    theme: cosmo
---

Working document. Tick items as they are completed; amend freely.

**Paper:** single combined manuscript — integrated population model for a breeding
Lesser Yellowlegs population, Anchorage, Alaska.
**Audience:** shorebird / ornithology researchers.
**Candidate journals:** *Ornithological Applications* (published the recent LEYE
Prairie Pothole migration paper, so the species and audience are established there),
*Ibis*, *Journal of Avian Biology*, *Wader Study*.

---

## 1. The argument

### Why anyone should care

- Lesser Yellowlegs declined ~63% over four decades and is a species of high
  conservation concern.
- **Field-based survival data for the species is scarce.** Existing population
  modeling rests on BBS/CBC indices fitted with state-space logistic models — not
  on measured vital rates.
- Harvest sustainability assessments in the Americas (18,000–47,000 birds/year)
  need adult survival estimates that currently barely exist.
- There is an active GPS-tracking programme for this species generating migration
  and connectivity papers. If tags materially reduce survival, that bears on those
  studies' inference and on deployment ethics.

### What this paper claims

1. Vital rates for a breeding LEYE population from an IPM combining counts,
   mark–resight, and productivity — adult apparent survival ≈ 0.63–0.70.
2. **Recruitment is delayed and low.** Locally-produced young recruit at ages 1 and 2
   in roughly equal proportion, not entirely at age 1 as commonly assumed. Now supported
   by **two independent eras** — 122 chicks banded 1995–99 and 133 banded 2018–25 —
   with near-identical age-at-first-return patterns thirty years apart. Modern total
   probability of recruiting locally ≈ 5%.
2b. **Immigration, not local production, appears to sustain the population.** With the
   modern recruitment prior in place, local recruits account for only ~32% of annual
   gains (down from 56% under the 1990s-derived prior), and immigration contributes
   ~15% of λ against recruitment's ~7%. This is the headline the tLTRE now delivers,
   and it has direct conservation implications: a local population propped up by
   immigration is vulnerable to declines elsewhere in the flyway.
   **State it with its range.** The ordering holds under every informative prior
   but the margin depends on which π the IPM receives: 15.4% vs 7.0% on the
   default, 13.1% vs 9.5% cohort-averaged, 15.5% vs 6.3% with the influential 1997
   cohort excluded — P(immigration > recruitment) 0.77, 0.62 and 0.84
   respectively. Claim the direction, not the ratio. See *Confounds to disclose*.
3. **GPS tags substantially reduce apparent survival** — a ~43% relative reduction,
   several times the published geolocator benchmark.
4. **Demographic stochasticity dominates variation in growth rate** in a population
   this small, which has direct consequences for monitoring design.
5. The population trend is negative but imprecise.

### Framing constraints — accept these early, do not overreach

- **λ = 0.930 (0.884–0.977), and the interval EXCLUDES 1.** P(λ < 1) = **0.997**,
  −7.0%/yr, and the fitted series falls to **0.436** of its 2014 level — a 56%
  decline over twelve years. **The decline can be claimed outright**, not hedged
  as "probably declining". Keep the hedges that are still real: abundance is
  relative, the recruitment/immigration split behind λ is not identified, and λ
  now describes the Anchorage bowl specifically rather than upper Cook Inlet.

  *Two revisions got this number here, and both belong in any methods narrative.*
  It read 0.959 (0.914–1.015), *including* 1 and "not a demonstrated decline",
  until **2026-11-03**, when that turned out to be the arithmetic mean of the
  annual rates mislabelled as a geometric mean — the geometric mean is the growth
  rate, and it is lower (0.939). On **2026-11-04** the count stream was restricted
  to the Anchorage bowl so that all three data streams describe the same
  population, which moved it again, to 0.930. Both changes made the decline
  steeper; neither was made in search of that. See *Confounds to disclose* for the
  arithmetic figure, which survives only as the tLTRE level decomposition's
  denominator, and *The count stream covers more ground* for the restriction.
- **Recruitment and immigration are confounded.** Both add birds and the relative-scale
  count index cannot separate them. Report the split as indicative.
- **Abundance is relative.** No absolute population size anywhere in the paper.
- **`psi` is not juvenile survival.** It is survive-*and*-return-locally. Do not
  compare it to published juvenile survival estimates; reviewers familiar with the
  literature will notice if we do.

### The tag result — how to make it credible

The GPS effect is large enough that reviewers will look for an artifact. Three
defences, in order of strength:

1. **Internal calibration.** Our geolocator estimate (≈10% absolute reduction)
   reproduces Weiser et al. (2016), who found ~8%/yr for leg-mounted geolocators
   across 16 Arctic-breeding shorebird species. The same model, same data, same
   priors returns the literature value for geolocators and a much larger value for
   GPS. That is hard to explain as a modeling artifact.
2. **Not a detection artifact.** The model estimates tag effects on resighting
   separately. GPS birds show no detection penalty (OR 0.92, P(<1) = 0.57), so the
   survival effect is not birds becoming harder to see. Notably geolocator birds
   *do* show one (OR 0.52, P(<1) = 0.97) — worth reporting in its own right.
3. **Mechanistically expected direction.** Brlík et al. (2020) found effects scale
   with relative device load and with harness attachment; GPS tags increase both.

**Caveat that must be stated plainly:** this is *apparent* survival. Tag-induced
dispersal is not separable from mortality here.

---

## 2. Status

> **Numbers in this document are a snapshot, and the fit moves.** Everything below was
> checked against `LEYE_IPM_samples.rds` on 2026-08-14. λ, the goodness-of-fit
> p-values and the sex effect all drifted slightly when the modern recruitment prior
> went in — small shifts, but the sex effect's lower bound moved onto exactly 1.00,
> which changes how it can be described. Re-verify every quoted value against the
> current posterior before drafting, and again before submission.
>
> **Values below were re-checked against `LEYE_IPM_samples.rds` on 2026-09-29**, after
> the refit on the real 1990s banding totals. Goodness-of-fit p-values are now
> counts 0.334, productivity 0.506, mark–resight 0.478.

| Component | State |
|---|---|
| IPM structure, fitting, convergence | Done |
| tLTRE | Done |
| Quarto report with publication-quality figures | Done |
| Goodness-of-fit / posterior predictive checks | Done — and it caught a real misspecification |
| Prior sensitivity | Done — and it **overturned** an earlier conclusion. Four-way (typical-cohort, cohort-averaged, 1997-excluded, vague); λ holds at **0.929–0.931** across all four, only the recruitment/immigration split moves. All four refitted 2026-11-04 on the restricted count stream |
| 2026 counts | **In** (added 2026-10-14). The series runs 2014–2026; the IPM fits 780 of the 935 records. λ 0.965 → 0.959 and P(λ<1) 0.88 → 0.94 — both **arithmetic-mean** figures on the full panel, as reported at the time; neither the pre-2026 nor the full-panel fit is the current model, so those values cannot be restated against it. The tLTRE split and every other conclusion held |
| λ summary | **Corrected 2026-11-03.** The reported growth rate was the arithmetic mean of the annual λ labelled as a geometric mean. Now the geometric mean everywhere it is a growth rate: on the full panel 0.959 → 0.939, P(λ<1) 0.94 → 0.997. The arithmetic mean survives only as the tLTRE level decomposition's denominator, where the additive identity requires it |
| **Spatial scope of the count stream** | **Changed 2026-11-04 — this is the current main model.** The counts are restricted to the 27 Anchorage-bowl wetlands so all three data streams describe one population. λ 0.939 → **0.930**, P(λ<1) = 0.997 either way, series ends at 0.436 of 2014 rather than 0.486. The wider-area fit is retained as the `full` sensitivity |
| Model structure in the report | Done — life-cycle figure, DAG and the full statistical formulation are now a section of the Quarto report, matching *Model at a glance* below |
| Tag-effect robustness analyses | **Re-run 2026-11-04** on the corrected 212-bird encounter history, max Rhat 1.036. GPS OR 0.28 (0.15–0.54), matching the IPM's 0.28 (0.15–0.53); duration and mass results unchanged to two decimals. One change: the transience offset's interval now **includes zero**. The morphometric check is now in the script rather than a console calculation, and reporting re-derives from the cache without a refit |
| Recruitment estimates | Now from the **modern era**; the 1990s cohorts still enter via the era comparison and shared detection — **confirmed 2026-11-04** that the 1990s chicks carried colour marks, which is what makes that sharing legitimate. Refitted 2026-09-29 on the real banding totals (122 chicks, 1995–1999) — no conclusion changed |
| Life-cycle figure and DAG | Done — both label the recruitment prior as the modern estimate, and name both cohorts (133 chicks 2018–25, 122 in 1995–99, with 178 1990s adults) |
| Novelty check (is this the first LEYE IPM?) | **Not done** |
| Study area map (F1) | **Done 2026-11-04** — `scripts/fig_study_area.R`, three panels (upper Cook Inlet, Anchorage bowl, Alaska locator) over a full-resolution GSHHG shoreline cached in `data/gis/`. **In the report** as `@fig-studyarea`, beside the data-sources table, using the two-panel `study_area_bowl` version (the 27 fitted wetlands plus an Alaska locator; the five unfitted wetlands are not drawn). The three-panel figure is still written for the manuscript. Still missing a mark–resight/nest-site layer: those coordinates are not in the repository. Drawing it **surfaced a confound nobody had recorded** — the count stream covered 67 km while the demographic data cover the Anchorage bowl; that is what led to the scope change above |
| Panel stability | **Tested 2026-11-04.** A nine-wetland balanced panel gives λ 0.950 with the interval crossing 1, on a third of the data and with a 52% wider interval — a bound on how much the claim leans on panel size, not a contradiction. `site_re` now monitored: it tracks observed richness (r = 0.98, slope 0.83), so the panel growing from 13 to 32 wetlands is not leaking abundance into the trend |
| Priors table (T1) | **Not done** |

---

## 3. Phase 0 — data resolution (resolved 2026-09-29)

**Resolved.** The real chick banding totals were recovered from the original records:
37 / 34 / 11 / 30 / 10 for 1995–1999, **122 chicks**, against the 139 the file used to
carry. 1995 and 1996 happened to be right; 1997 and 1998 had been simulated and were
overstated by 9 and 18, and the 1999 cohort was missing altogether. The whole chain was
rebuilt and refitted on that day.

- [x] Recover the real chick banding totals from the original records.
- [x] Rebuild `data/LEYE_juv_EH_1995_2000.csv` from real denominators. The simulation
      branch is deleted; cohort sizes are a documented constant with a `stopifnot`
      guard, so the file no longer depends on the seed.
- [x] Refit `scripts/LEYE_recruitment_two_era.R` → new `pi_M` / `kappa` priors.
- [x] Update the Beta priors in `scripts/LEYE_IPM.R`
      (`dbeta(1.8410, 32.5058)` and `dbeta(7.6338, 8.6569)`, from the modern era).
- [x] Refit the IPM, re-run the tLTRE, re-render the report.
- [x] **Re-check whether any conclusion changed.** None did — see below.
- [x] **Confirm whether the 1990s chicks carried colour marks.** Confirmed 2026-11-04
      (J. Wagner, from the field records). The juvenile module fits one detection
      process across 1990s chicks and the 178 contemporaneous adults, and that is
      only legitimate if both were read the same way. The justification had been
      circumstantial — nicknames like "Firehouse Blue" look like colour-mark
      combinations — and the alternative, metal bands only, would have meant lower
      true chick detection and an **underestimated** `psi`. Nothing is refitted: the
      module already assumed shared detection, so the estimates do not move. What
      changes is that the assumption is now supported rather than inferred, which is
      the difference between a caveat and a methods sentence.

**What moved.** Smaller 1990s denominators raise π₉₀ from 0.115 to 0.139, and because
detection is shared across eras that carries into the modern estimate too: π_M 0.049 →
0.054, κ 0.485 → 0.469. Downstream the changes are small — λ 0.963 → 0.965 (arithmetic
mean, the summary in use at the time; see the λ-summary correction below), the tLTRE
level split 76/17/7 → 75/17/7, P(immigration > recruitment) 0.83 → 0.81, GPS odds ratio
0.29 → 0.28. **Every conclusion in this document survives.** That robustness is itself
worth a sentence in the manuscript: the headline did not depend on the simulated
denominators, which is not something we could assert before.

*(Those are the values as they stood on 2026-09-29, computed on the full 32-wetland
count stream and with the arithmetic λ. They are kept as the record of what that change
did; for current values see the status table. On the present model the split is
77/15/7 and P(immigration > recruitment) is 0.77.)*

**The one thing to disclose.** The 1997 cohort is 11 chicks with 7 returning — a 64%
local return rate against 11%, 7%, 0% and 0% in the other 1990s cohorts and 4% in the
modern era. 1997 was a genuinely low-effort season (21 adults marked, against
52/48/36/21 around it), so the small denominator is real, but one cohort now carries
most of the 1990s recruitment signal. Say so rather than leaving a reviewer to find it.

**Still blocking the era comparison**, though for a different reason than before: the
contrast is now P(decline) = 0.91 with a ratio spanning 0.08–1.59. That is suggestive,
not a result, and the width comes from 13 and 5 returning birds — not from the
denominators, which are now real.

The contrast is now drawn as **F7** and reproduced in the report under
*Has recruitment changed between eras?*, built from the module posterior at
render time rather than embedded as an image. The figure is deliberately
two-panel: the era estimates alone read as a clean halving, and the ratio
posterior beside them is what stops that reading. Use both or neither.

---

## 4. Phase 1 — analyses to add before writing

### Reviewers will ask for these

- [x] **Posterior predictive checks / goodness-of-fit.** Done — Freeman-Tukey
      discrepancies per data stream, built into `LEYE_IPM.R` section 2f.
      **This found a real misspecification.** The original Poisson on chicks-per-nest
      failed outright (Bayesian p = 0.000): the data are bimodal (0, 2, 3, 4 chicks;
      *no* nests with exactly 1) because a 4-egg-clutch shorebird either loses the
      nest or fledges a near-full brood, and a Poisson cannot produce that shape.
      Note `var/mean = 1.18` — a routine overdispersion check would have passed it.
      Replaced with nest success × Binomial(4, q). Post-fix p-values:
      **counts 0.313, productivity 0.498, mark-resight 0.475.**
- [x] **Prior sensitivity on recruitment.** Done, via prior hyperparameters passed
      as constants (so each sensitivity run is a config change, not a forked model).
      **This conclusion REVERSED when the modern recruitment prior replaced the
      1990s one, and the earlier version of this entry was wrong.**
      Now a **four-way** comparison — typical-cohort (default), cohort-averaged,
      1997-excluded, and vague — reported in the Quarto report's prior-sensitivity
      section with π, κ, ω, λ and the contribution shares for each. All three
      companion fits are built by `scripts/LEYE_IPM_sensitivity.R` (no argument
      runs all three; `cohortmean`, `drop97` or `vague` runs one; `restamp`
      rewrites provenance without refitting), which takes the data, structure and
      sampler settings from `scripts/LEYE_IPM.R` and overrides only the
      recruitment hyperparameters — so the companions stay comparable to the main
      fit by construction rather than by anyone remembering to keep them in step.
      The vague fit previously had no script at all; it was produced by
      hand-editing constants, with no record of how. Each is checked against the
      main fit's *structural* fingerprint at render time, so a stale companion is
      flagged in the report rather than silently compared.

      Under the 1990s prior, informative and vague agreed (0.106 vs 0.108) and the
      recruitment level looked safely data-driven. Under the modern prior they
      differ two-fold — `pi_rec` = **0.052 informative vs 0.105 vague** — because
      the counts and the modern chick resights disagree, and the tighter modern
      prior wins. The old agreement was a coincidence: the 1990s prior happened to
      sit near what the counts wanted.

      The **1997-excluded** setting (added 2026-11-03) probes the one direction the
      other two did not. Both of those shift π's centre and leave its width alone;
      this one makes the prior *tighter* (sd 0.021 vs 0.038), because dropping the
      1997 cohort collapses σ_π from 1.39 to 0.44 and with it the cohort-level
      heterogeneity. Both recruitment parameters are re-derived from the same
      re-fit, so κ moves too (0.47 → 0.40) rather than mixing eras of evidence. It
      is the **most adverse** of the four for the headline: π 0.045, ω 0.143,
      recruitment 6.3% of λ against immigration's 15.5%, and
      P(immigration > recruitment) = 0.84 — the ordering gets *stronger*, not
      weaker, when the influential cohort is removed. λ is 0.929.

      **This belongs in the Discussion, not buried in a supplement.** Three readings
      cannot be separated with current data: (a) the counts cannot distinguish
      recruitment from immigration, so count-inferred recruitment absorbs
      unmodeled inflow; (b) assuming modern chick detection matches the 1990s may
      be too generous, biasing recruitment down; (c) recruitment genuinely is low
      and immigration genuinely high.

      **Timing** is still prior-driven: under vague priors `kappa`'s interval
      widens (0.24–0.70) to (0.02–0.97) with its centre unmoved, so the counts say
      nothing about the age split. Its centre moves only in the 1997-excluded fit
      (0.465 → 0.401), which is the only setting that also replaces κ's prior — the
      chick data, not the counts, are doing that.

      **Convergence, checked and fixed 2026-11-03.** The recruitment parameters
      are clean in all four fits — Rhat ≤ 1.003 for π and κ, ≤ 1.010 for ω. The
      earlier "Rhat ≤ 1.004" claim in this plan covered only those and should not
      be quoted as a whole-model figure, because `logB0` reached **1.188 in the
      vague fit**. The vague fit was therefore re-run at 300k iterations (3×,
      59 min), which fixed it: `logB0` Rhat 1.188 → **1.008**, ESS 110 → 325, and
      the posterior did not move (π 0.103 → 0.105, κ 0.484 → 0.484, ω 0.114 →
      0.113, λ unchanged to four decimals). `scripts/LEYE_IPM_sensitivity.R` now
      takes `LEYE_NITER` / `LEYE_NBURNIN` / `LEYE_THIN` from the environment so a
      longer chain needs no forked script, and the report now prints each
      companion's worst Rhat beside its priors, since the convergence table covers
      the main fit only and a badly mixed companion was otherwise tabulated with
      no warning.

      **`logB0` is weakly identified in every fit, by construction** — the count
      model is `log(λ_i) = logB0 + log(N_tot) + …`, so only the sum is pinned, and
      cor(logB0, mean log N_ad) = −0.85 everywhere. Its ESS is ~90–150 out of
      21 000 draws in every fit, so at that ESS its Rhat is itself a noisy
      statistic and lands anywhere from 1.01 to 1.19 on luck alone. **The current
      main fit sits at logB0 Rhat 1.065 with ESS 153** — a higher Rhat than the
      previous main fit (1.020) on *better* mixing (ESS 89), which is the clearest
      demonstration available that the statistic is noise here. Everything else in
      that fit is clean: excluding `logB0`, max Rhat ≤ 1.026 over λ, N_tot, N_ad
      and the vital rates. It is not a vague-prior pathology, and the
      report's *Convergence* section already explains the structural reason. A
      slice sampler on `logB0` was tested and rejected: ESS per 1000 iterations
      goes 0.5 → 1.4, a 3× gain where ~100× is needed, and nothing for π or ω
      (ESS 303 → 275, 255 → 245). No univariate sampler can traverse this ridge,
      because moving along it requires `logB0` and the whole discrete `N_ad`
      trajectory to shift together. **Longer chains are the only lever** — ESS
      scales linearly with iterations, as the 3× run confirmed.
- [ ] **Novelty check.** Establish whether a field-based IPM exists for this species.
      Frame claims accordingly.

### Strengthening the tag result

**Every number in this subsection was re-derived on 2026-11-04** from the corrected
212-bird encounter history (`Rscript scripts/LEYE_tag_effects_CJS.R`, three models,
max Rhat 1.036 on the quantities quoted). The previous values came from a 2026-08-06
fit of the 213-bird history and are given below only where they moved.

- [x] **Time-since-tagging effect.** Done — `scripts/LEYE_tag_effects_CJS.R`.
      Duration slope −0.091 (−0.401, 0.219), P(<0) = 0.72: direction matches
      Pakanen but no detectable effect, and WAIC worsens (848.7, against 845.6 for
      M0 and 844.8 for M1). The GPS effect survives its inclusion — 0.29
      (0.12–0.65), P(<1) = 0.998.
- [x] **Selection bias.** Done. **GPS birds were significantly heavier at banding:
      87.3 g vs 82.0 g untagged (+5.31, 95% CI 3.08–7.54, p < 0.0001).** Geolocators
      were not (−0.22, p = 0.84). This works *against* the finding — heavier birds
      should survive better — so the tag effect is conservative. This comparison used
      to be a console calculation with no script, which is why it could not be
      re-derived when the encounter history changed; it is now Section 1B of
      `LEYE_tag_effects_CJS.R` and runs off the same objects the models are fitted to.
- [x] **Standalone CJS check.** Done, and **both sides are now current**. M0 gives GPS
      OR 0.28 (0.15–0.54) against the IPM's 0.28 (0.15–0.53) — the tag effect is not
      an artifact of IPM linkage. The geolocator estimates agree less closely, 0.53
      (0.30–0.92) standalone against 0.65 (0.37–1.13) in the IPM, but the intervals
      overlap across almost their whole length and both sit below 1.
- [x] **Transience control.** Added because exposure differs by tag class: untagged
      birds carry 16.6% of their monitored intervals as first intervals, against 8.4%
      for geolocator and 13.8% for GPS birds. Including it does not weaken the GPS
      effect (0.28 → 0.26). **Note the denominator** — see the exposure bullet under
      *Confounds* below, because the ordering reverses under the other obvious one.
- [ ] **Tag mass as % of body mass.** Use **87.3 g** as the GPS denominator, not the
      82.7 g population mean — GPS birds were a heavier subset. Still need tag specs.
- [x] **Re-run the tag analyses against the corrected encounter history.** Done
      2026-11-04. The change was one bird out of 213, neither one tagged, and the
      estimates duly barely moved: GPS OR 0.28 in M0 either way, duration slope
      −0.090 → −0.091, mass difference +5.33 → +5.31 g. **One thing did move enough
      to change a sentence**: the transience offset's 95% interval now includes zero
      (below). Reporting is now re-derivable from the cache without a refit —
      `SKIP_TAG_RUN <- TRUE; source("scripts/LEYE_tag_effects_CJS.R")` reprints every
      number in seconds, which is what the three-month-stale set was missing.
- [ ] Attachment method (harness vs leg-mount) for each tag type — a documented
      effect modifier (Brlík et al. 2020).
- [ ] **Tag recovery dates**, if they exist. Geolocators must be recovered, so some
      birds stopped carrying tags mid-study; the duration test currently uses
      time-since-banding as an attenuated proxy. Recovery dates would sharpen it.

### Confounds to disclose in the manuscript

Found during the robustness work; all belong in the Methods or Discussion.

- **The three data streams do not cover the same ground.** Counts span 67 km and 32
  wetlands; survival and productivity come from the Anchorage bowl. Tested rather than
  assumed — restricting the counts to the demographic footprint makes the decline
  *steeper*, not weaker. Full numbers in *The count stream covers more ground than the
  demographic data*, below.
- **Tag type is largely confounded with banding era** — geolocators 2010–2018,
  GPS 2018–2020, untagged 2011–2025 in the analysis window. Not fatal: 52 untagged
  birds were banded in the same 2018–2020 window as the GPS birds, so the contrast
  is not purely across eras. Year-varying baseline survival absorbs period effects.
- **Sex composition differs** — GPS 58.7% female against 41.0% of untagged birds,
  and no GPS bird is of unknown sex while 18.1% of untagged birds are. Since males
  survive better (OR 1.51, 0.99–2.32), this biases GPS *downward*; it is controlled
  by the sex covariate, but the model must retain sex for that reason. *(The 44.8%
  previously quoted here is not reproducible from the analysis set under any
  definition tried; 41.0% is what Section 1B computes from the birds the models are
  fitted to.)*
- **Transience runs the other way than expected, but it is no longer a result** —
  the first-interval offset is *positive*, +0.638 (−0.027, 1.417), i.e. survival
  looks higher just after banding, with P(> 0) = 0.97. **On the corrected encounter
  history the 95% interval includes zero**, where the pre-correction fit gave
  +0.659 (0.003, 1.417) and excluded it. One bird moved it across the line, which
  is the honest description of how thin the evidence was all along. Report it as a
  direction worth a sentence, not an effect, and do not call it significant.
- **Exposure imbalance depends on the denominator, and only one of them is the
  right one.** Counting every interval the likelihood uses — banding to the end of
  the study, which is what the marginalised CJS integrates over — untagged birds
  carry 16.6% of their intervals as first intervals, against 8.4% for geolocator
  and 13.8% for GPS birds, so the untagged class is the one most exposed to a
  first-interval effect and M1 is the control for it. Counting only intervals up to
  each bird's last detection reverses the order completely (untagged 35.6%, geo
  48.8%, GPS 64.8%), because GPS birds are resighted least. The transience argument
  is only valid under the first denominator; it is pinned in code in Section 1B of
  `LEYE_tag_effects_CJS.R` so the manuscript cannot quietly pick up the other one.
- **Recruitment is prior-driven, and immigration absorbs the difference**
  (re-checked 2026-11-04 on the current fit). Under the default prior the
  posterior barely leaves it: π moves 0.0536 → **0.0508**, a shift of −0.0028 or
  **−0.07 prior SD**, and the data remove only ~10% of the prior spread
  (interval 0.006–0.149 → 0.006–0.134). For κ there is no information at all:
  a shift of −0.02 prior SD and an interval that stays the same width to three
  decimals (0.464 → 0.466). Under vague priors κ returns 0.484 (0.023–0.971),
  which *is* Beta(1,1).

  **The report said the opposite until 2026-11-04.** The callout under the
  prior-versus-posterior figure claimed `pi_rec` "moves appreciably from its
  prior", which was never true of any fit in this repository — it was the kind
  of statement that sounds like it was checked. It and the README's copy are
  corrected, and the report now computes the shift and the narrowing inline from
  the fit rather than describing them in words, so the claim cannot drift from
  the model again. The *direction* of the original point survives — the counts
  say more about the level than the timing — but "appreciably" was wrong, and
  the honest summary is that they say very little about either.

  λ below is the **geometric** mean across years, the growth rate (see the next
  bullet). The shares are per-capita contributions with their own denominator,
  so they are unaffected by that choice.

  All four refitted 2026-11-04 on the restricted (27-wetland) count stream, so the
  columns are comparable to the main model and to each other.

  | prior on π | π | κ | ω | λ | rec. share | imm. share | P(imm > rec) |
  |---|---|---|---|---|---|---|---|
  | typical cohort (default) | 0.051 | 0.466 | 0.143 | 0.930 | 7.0% | 15.4% | 0.77 |
  | cohort-averaged | 0.069 | 0.465 | 0.121 | 0.931 | 9.5% | 13.1% | 0.62 |
  | 1997 cohort excluded | 0.045 | 0.400 | 0.143 | 0.929 | 6.3% | 15.5% | 0.84 |
  | vague | 0.096 | 0.484 | 0.108 | 0.931 | 12.8% | 11.5% | 0.46 |

  The four runs trace the recruitment–immigration ridge: as π rises, ω falls to
  compensate and **λ barely moves** — 0.929–0.931, a spread of 0.002 against a
  95% CRI nearly fifty times that wide (≈ 0.88–0.98). The paper's growth-rate
  results are immune to this choice; only the *attribution* is not.
  Disclose it rather than let a reviewer find it: it is the empirical signature
  of the additive ψ/ω confounding set out in *Model at a glance*, and it is
  precisely why π carries an informative prior. The honest framing is that the
  chick data identify recruitment and the counts identify λ; neither alone
  separates recruitment from immigration.

- **Two averages of λ exist in this analysis, and the paper reported the wrong
  one as the growth rate** (found and fixed 2026-11-03). The reported figure was
  the **arithmetic** mean of the twelve annual λ, averaged across draws, but it
  was labelled a geometric mean in the report's KPI card, in the λ figure's
  subtitle, in the prior-sensitivity table and in this plan.

  Values below are from the current (27-wetland) main fit.

  | summary | value | 95% CRI | P(λ<1) | implied 12-yr change |
  |---|---|---|---|---|
  | geometric mean (correct) | 0.930 | 0.884–0.977 | 0.997 | −58% |
  | arithmetic mean (was reported) | 0.953 | 0.903–1.018 | 0.936 | −44% |

  The geometric mean is provably the right one: within each posterior draw it
  equals (N_tot[2026]/N_tot[2014])^(1/12) to machine precision — max absolute
  difference 2.2 × 10⁻¹⁶ across 21 000 draws, correlation 1.000000, against 0.898
  for the arithmetic version. It is an algebraic identity, not an approximation.
  The fitted series ends at 0.436 of its 2014 level, i.e. −56.4%, which the
  geometric figure reproduces and the arithmetic one does not. (The −58% in the
  table is mean(λ)¹² and the −56.4% is mean(λ¹²); they differ for the same
  reason the two means do, and −56.4% is the one to quote.)

  **Why the two differ, in one paragraph** (this is the "Jensen" referred to
  throughout this document). Growth rates multiply rather than add, so averaging
  them arithmetically is the wrong operation. Gain 20% one year and lose 20% the
  next and the two λ, 1.2 and 0.8, average to exactly 1.0 — but the population
  ends at 1.2 × 0.8 = 0.96 of where it started, a true rate of √0.96 = 0.98 a
  year. The arithmetic mean of a set of ratios is never below their geometric
  mean, and the gap widens with the spread: to second order about σ²/2 on the
  log scale, which here predicts 0.0229 against the 0.0229 observed (the mean
  within-draw SD of log λ is 0.221). **Jensen's
  inequality** is the general statement — put a quantity through a curved
  transformation (here a logarithm) and averaging before is not the same as
  averaging after, with the direction set by which way the curve bends. The
  same mechanism, with ilogit as the curve instead of log, is what separates the
  typical cohort's π from the cohort average two bullets down.

  **This makes the decline claim stronger, not weaker** — the 95% CRI moves from
  straddling 1 to lying wholly below it — so it is a correction to make before
  submission, not a problem to manage. The *Framing constraints* bullet that told
  us to hedge has been rewritten accordingly.

  The arithmetic mean is still correct in exactly one place: the tLTRE **level**
  decomposition, where λ_s = survival + recruitment + immigration holds exactly
  each year and only the arithmetic mean preserves that identity across years, so
  the pathway shares sum to 100%. That figure also covers 11 of the 12
  transitions, because the decomposition needs the previous year's productivity,
  and the transition it drops (2014→2015, λ = 1.139) is the second highest in
  the series — so dropping it lowers the window's mean. Hence three different
  numbers, all correct for their own purpose: 0.930 (growth), 0.953
  (whole-series arithmetic), 0.936 (tLTRE window arithmetic, the level figure's
  denominator; geometric on the same window is 0.916). The report spells this
  out in a callout under that figure, and
  `tltre_level()` carries a comment telling the next person not to "fix" it.

- **Which π the IPM should receive is a genuine modeling choice** (added
  2026-10-02). `LEYE_recruitment_two_era.R` fits a cohort random effect on the
  logit scale and reports π_M = ilogit(μ_π + b_era) — the *typical* cohort, with
  its random effect at zero. But the IPM uses π in an expected-count expression,
  R = N·(f/2)·ψ, and has no cohort random effect of its own, so arguably it wants
  E[π_c], the mean **across** cohorts. With σ_π ≈ 1.39 the two differ by about a
  third (0.054 vs 0.072), because ε_c is symmetric about zero on the *logit*
  scale but ilogit bends: cohorts above the typical one land further from it in
  probability than cohorts below, so back-transforming the typical cohort gives
  a smaller number than averaging the back-transformed cohorts. Same mechanism
  as the λ averaging above (Jensen's inequality), with ilogit as the curve
  instead of log. Both are now fitted and
  reported side by side; neither is obviously wrong. Hold the cohort-averaged
  version loosely, because σ_π rests on 13 cohorts and 18 returning birds and the
  1997 cohort alone sits at π_c = 0.58 against 0.03–0.18 elsewhere — and since
  σ_π is shared across eras, that one cohort inflates the gap for the modern
  cohorts too. **Report the immigration ordering at the conservative end:**
  immigration exceeds local recruitment with P ≈ 0.67–0.88 depending on the
  definition, a margin of roughly 1.5:1 to 2.8:1 rather than a settled number.

  That σ_π caveat is now **quantified** (2026-11-03).
  `scripts/LEYE_recruitment_drop_cohort.R 1997` refits the two-era module with the
  1997 cohort removed, inheriting every other prior, the detection process, adult
  survival, κ and the likelihoods unchanged. One cohort is carrying most of the
  era contrast: π_90 0.139 → 0.074, b_era −1.125 → −0.560, P(declined) 0.910 →
  0.809. The decisive number is σ_π, which falls from 1.39 to 0.44 —
  P(σ_π < 0.5) goes from 0.012 to **0.658**. Without 1997 there is no detectable
  cohort-level heterogeneity anywhere in the dataset, which also means the
  cohort-averaged and typical-cohort definitions of π nearly coincide and the
  typical-vs-average gap that motivates the cohort-averaged setting largely
  disappears — with σ_π near zero there is no spread for the curve to act on. Report
  the era decline with the leave-one-out value alongside it; it is not a result
  that should be stated as 0.91 alone.

### Spatial scope of the count stream — restricted to the Anchorage bowl (2026-11-04)

**Surfaced by drawing F1, and it changed the main model.** The three data streams did
not describe the same place. Every marked bird and every monitored nest comes from the
Anchorage bowl. The counts came from 32 wetlands spanning **67 km**, and five of them —
both Palmer Hay Flats accesses, Eklutna, Big Lake and Lynn Fuller's Ponds, 37–67 km
out — supply **22.6% of all birds counted** on 16.6% of the surveys. Nothing in the
repository said so before F1 was drawn.

**The IPM now fits the 27 bowl wetlands** (`IPM_COUNT_WETLANDS` in `LEYE_IPM.R`, 780 of
935 records), so λ and the vital rates that are supposed to explain it describe one
population. The wider-area fit is retained as the `full` sensitivity.

A second problem sits underneath it. The panel is not stable: it grew from **13
wetlands in 2014 to 32 by 2021**, only **9** were surveyed in all 13 years, and the
wetlands added late are the rich ones — correlation between a site's mean count and
its mean survey year **+0.50**, with the two richest of all (Coastal Refuge 7.8
birds/survey, Oceanview Ghost Forest 5.2) surveyed in only 7 and 8 years centred on
2023.

The model is not naive about this: every wetland already gets its own level offset,
`site_re`, with σ_site ≈ 1.5 on the log scale. What it *assumes* is that the wetlands
it keeps share one temporal trajectory, `N_tot[t]`. Three fits, differing only in which
wetlands enter the count likelihood — same priors, chain length, seeds and monitors:

| fit | count stream | records | geometric-mean λ | P(λ<1) | N₂₀₂₆/N₂₀₁₄ |
|---|---|---|---|---|---|
| **main** | Anchorage bowl, 27 wetlands | 780 | **0.930** (0.884–0.977) | 0.997 | 0.436 |
| `full` | all 32 wetlands | 935 | 0.939 (0.897–0.982) | 0.997 | 0.486 |
| `balanced9` | surveyed in all 13 years | 333 | 0.950 (0.886–1.015) | 0.937 | 0.583 |

**Restricting the scope was not a way of removing an inconvenient signal — it makes
the decline steeper.** λ 0.939 → 0.930, with P(λ<1) = 0.997 either way. The distant
wetlands were damping the decline slightly, not creating it. Say this explicitly in
the manuscript: a reviewer's first thought on seeing sites dropped is that dropping
them helped the result, and here it did the opposite. It also means there is no
divergence to model, so the group × year contrast sketched as a contingency should
**not** be added — the test for it came back negative.

The change was made because the alternative is a λ that describes a different
population from the rates decomposing it, which would make the tLTRE meaningless, not
because of what it does to the number.

**`balanced9` is the honest limit on the claim.** It flattens λ to 0.950 and its
interval crosses 1 (P = 0.937). But it is a third of the data and its interval is 52%
wider (0.129 vs 0.085), so the shift in the point estimate sits well inside the
uncertainty of either fit; σ_site also inflates from 1.52 to 2.21 on nine wetlands.
This is not evidence against the decline. It *is* the reason not to write "the
interval excludes 1" as though it were assumption-free: that holds on the main fit and
on the full panel, not on the nine-wetland panel.

**The demographic rates barely notice.** π 0.051/0.051/0.052, κ 0.466/0.465/0.467,
φ_ad 0.699/0.701/0.701 across main/`full`/`balanced9`; only ω moves at all
(0.143/0.157/0.156). The count stream is shaping the latent trajectory, not pulling
the vital rates around — which is the reassurance that restricting it has not quietly
changed what the mark–resight and nest data are saying.

**`site_re` was monitored for this and answers the shrinkage worry directly.** Across
the 29 wetlands with any detection, the fitted site effect tracks observed richness
almost exactly: **cor(site_re, log mean count) = 0.981, slope 0.83**. The 17%
compression is ordinary partial pooling, nowhere near enough to push a rich site's
abundance into `N_tot` for the late years it was surveyed. An independent check agrees:
`balanced9` *excludes* Coastal Refuge, the one thin-rich-late site that could plausibly
leak, and λ gets flatter rather than steeper — the opposite of what the leak would
predict.

**One prediction here was wrong and is worth recording.** Crude log-linear year slopes
(−0.011 raw, −0.019 balanced, −0.029 with site fixed effects) suggested the balanced
panel would come out *steeper*. It came out flatter. The crude slopes had no site
effects; the IPM already has them, so there was nothing for a balanced panel to repair.
Do not quote those slopes — they describe a model this project does not fit.

**Convergence.** Excluding `logB0`, max Rhat ≤ 1.026 in all three. The main fit shows
logB0 Rhat 1.065, but its ESS is 153 against `full`'s 89 — that is the documented
weak-identifiability ridge producing a noisy statistic, not a consequence of the
restriction.

**For the manuscript:** state in Methods that the count stream is restricted to the
Anchorage bowl and why, report `full` as the wider-area sensitivity with the direction
of the change, and say plainly that λ, survival and productivity now all describe the
same ~25 km area. What remains unprovable either way is connectivity: the encounter
histories cannot show movement between the bowl and Palmer Hay Flats, because
`Study_Region` records where a bird was *banded*, not where it was resighted. So the
restriction cannot be defended as "those are a separate population" — only as "those
are not the population we have demographic data for", which is the honest claim and
the sufficient one.

Scripts: `LEYE_IPM_sensitivity.R full balanced9`, compared against the main fit by
`scripts/compare_count_subsets.R` → `output/plots/count_subset_comparison.png`.

### Optional but valuable

- [x] Balanced-panel count sensitivity — **done 2026-11-04 inside the IPM**, as the
      `balanced9` setting above, which is the stronger version of what
      `scripts/LEYE_trend_NB_coreSites_balancedPanel.R` does outside it. Coverage grew
      13 → 32 wetlands; see the subsection above for what that does and does not do
      to λ.
- [ ] Power/monitoring framing: given demographic stochasticity dominates Var(λ),
      how long must this survey run to detect a decline of a given size? This turns
      a limitation into a contribution for a monitoring-oriented audience.

---

## 5. Phase 2 — structure

### Methods

1. **Study area and field methods** — Anchorage wetlands; volunteer survey design
   (3 weekly visits, 32 wetlands surveyed, the 27 in the Anchorage bowl used by the
   IPM); banding and resighting; nest monitoring; tag
   deployment (types, masses, attachment, years).
2. **Data sources** — counts 2014–2026 (935 records, 780 fitted); adult encounter histories
   2010–2026 (**212** individuals); productivity 2018–2026 (77 nests); chick
   cohorts from **two eras** — 122 banded 1995–99 and 133 banded 2018–25.

   The adult total fell 213 → 212 because band 129232054, banded as a chick in
   2023 and recaptured as a breeding adult in 2025, was being counted twice: once
   as a chick and again as a new adult banding. It is a confirmed local recruit
   and now sits only in the chick history. Worth a sentence in Methods — it is the
   only directly observed recruitment event in the modern data.
3. **IPM structure** — process model with delayed recruitment; the two-year lag and
   why; survival with sex and tag covariates; relative-scale count observation model.
4. **Recruitment sub-model** — the two-era chick+adult joint analysis, marginalising
   over recruitment age. Explain *why*: a chick first seen at age 2 may have recruited
   at age 1 and been missed, so raw return rates confound timing with detection.
5. **tLTRE** — level and variance decomposition; explicit separation of demographic
   stochasticity.
6. **Implementation** — NIMBLE, chains/iterations/convergence, priors table.

---

### Model at a glance

Transcribed from `scripts/LEYE_IPM.R`, re-checked 2026-10-14. Indices: *t* = year (1…T, 2014–2026),
*i* = count record, *j* = marked individual, *m* = nest, *w* = survey week, *k* = wetland.
**If the model changes, this block must be re-transcribed — it is documentation, not code.**

**Process model.** Two free initial years, because recruitment lags two years and
there is nothing to build $N_1$ and $N_2$ from:

$$
\begin{aligned}
N_1,\, N_2 &\sim \mathrm{Normal}(0,\ 50)\ \text{truncated at } 0 \\[2pt]
R_t &\sim \mathrm{Poisson}\!\left(N_{t-1}\,\frac{f_{t-1}}{2}\,\psi_1
        \;+\; N_{t-2}\,\frac{f_{t-2}}{2}\,\psi_2\right) \\[2pt]
S_t &\sim \mathrm{Poisson}\!\left(N_{t-1}\,\phi_{t-1}\right) \\[2pt]
I_t &\sim \mathrm{Poisson}\!\left(N_{t-1}\,\omega\right) \\[2pt]
N_t &= S_t + R_t + I_t \\[2pt]
\lambda_t &= N_{t+1} \,/\, N_t
\end{aligned}
\qquad t = 3,\dots,T
$$

Recruits enter the breeding class directly; there is no juvenile stage.

**Vital rates.**

$$
\begin{aligned}
\operatorname{logit} \phi_t &= \mu_\phi + \varepsilon_{\phi,t},
  &\varepsilon_{\phi,t} &\sim \mathrm{Normal}(0,\ \sigma_\phi)
  &&\text{adult survival} \\[2pt]
\operatorname{logit} s_t &= \mu_s + \varepsilon_{s,t},
  &\varepsilon_{s,t} &\sim \mathrm{Normal}(0,\ \sigma_s)
  &&\text{nest success} \\[2pt]
f_t &= s_t \cdot C \cdot q &&&&\text{fledglings per breeding pair} \\[2pt]
\psi_1 &= \pi\,\kappa &&&&\text{recruit at age 1} \\[2pt]
\psi_2 &= \pi\,(1-\kappa) &&&&\text{recruit at age 2}
\end{aligned}
$$

with $C = 4$ (clutch size) and $q$ the per-egg fledging probability. The
immigration rate $\omega$ is constant — see the note at its prior.

**Observation 1 · counts** (relative index):

$$
\begin{aligned}
\log \mu_i &= \beta_0 + \log N_{t(i)} + \text{week}_{w(i)}
   + b_e E^{(1)}_i + b_d E^{(2)}_i + u_{k(i)} + \delta_i \\[2pt]
\text{count}_i &\sim \mathrm{Poisson}(\mu_i) \\[2pt]
u_k &\sim \mathrm{Normal}(0,\ \sigma_{\text{site}})
   &&\text{wetland random effect} \\[2pt]
\delta_i &\sim \mathrm{Normal}(0,\ \sigma_{\text{od}})
   &&\text{overdispersion}
\end{aligned}
$$

$\beta_0$ is non-identifiable against the level of $N$ — read shape, not magnitude.

**Observation 2 · mark–resight** (CJS with latent alive-state $z$):

$$
\begin{aligned}
\operatorname{logit} \phi_{j,t} &= \ell_\phi\!\left[\mathrm{lookup}(t)\right]
   + \beta^{\text{sex}}_{\mathrm{sex}(j)} + \beta^{\text{tag}}_{\mathrm{tag}(j)} \\[2pt]
\operatorname{logit} p_j &= p_0 + \beta^{p}_{\mathrm{tag}(j)} \\[2pt]
z_{j,t} &\sim \mathrm{Bernoulli}\!\left(\phi_{j,t-1}\, z_{j,t-1}\right),
   \qquad z_{j,\,f(j)} = 1 \\[2pt]
y_{j,t} &\sim \mathrm{Bernoulli}\!\left(p_j\, z_{j,t}\right)
\end{aligned}
$$

$$
\ell_\phi[1{:}T-1] = \operatorname{logit}\phi_t, \qquad
\ell_\phi[T] = \operatorname{logit}\phi_{\text{pre}}, \qquad
\ell_\phi[T+1] = \operatorname{logit}\phi_{\text{post}}
$$

The two overflow slots absorb intervals falling outside the count series.

**Observation 3 · productivity** (two-part):

$$
\begin{aligned}
\text{succ}_m &\sim \mathrm{Bernoulli}\!\left(s_{t(m)}\right)
   &&\text{did the nest fledge anything?} \\[2pt]
\text{chicks}_m &\sim \mathrm{Binomial}\!\left(C,\ q \cdot \text{succ}_m\right)
   &&\text{how many of the clutch}
\end{aligned}
$$

A single Poisson was **rejected** here ($p = 0.000$): the data are bimodal
(0, 2, 3, 4 chicks, none with exactly 1).

**Priors.**

$$
\begin{aligned}
\mu_\phi,\ \mu_s,\ p_0,\ \beta^{\text{sex}},\ \beta^{\text{tag}},\ \beta^{p}
   &\sim \mathrm{Normal}(0,\ 1.5) \\[2pt]
\beta_0,\ b_e,\ b_d,\ \text{week}_w &\sim \mathrm{Normal}(0,\ 5) \\[2pt]
\sigma_\phi &\sim \mathrm{Exponential}(2) \\[2pt]
\sigma_s,\ \sigma_{\text{site}},\ \sigma_{\text{od}},\ \omega
   &\sim \mathrm{Exponential}(1) \\[2pt]
\phi_{\text{pre}},\ \phi_{\text{post}},\ q &\sim \mathrm{Uniform}(0,\ 1) \\[2pt]
\pi &\sim \mathrm{Beta}(1.8410,\ 32.5058)
   &&\text{modern-era recruitment, external module} \\[2pt]
\kappa &\sim \mathrm{Beta}(7.6338,\ 8.6569)
   &&\text{age-1 fraction, external module}
\end{aligned}
$$

$\pi$ and $\kappa$ are the **only** informative priors, and both are passed as
constants so the vague-prior sensitivity run is a configuration change rather
than a model edit.

**Reading the structure.** Three things are unusual and each should be defended
explicitly in Methods rather than left for a reviewer to find:

- **Nₜ is breeding adults only.** Recruits join it directly and there is no observable
  juvenile class, because the counts cannot separate non-breeding age-1 birds from
  adults. This is why ψ is a *composite* survive-and-return probability, not survival.
- **ψ and ω enter λ additively** (`λ ≈ φ + ω + ψ₁f/2 + …`), so only their sum is well
  determined by a relative-scale index. That is the source of the recruitment/
  immigration confounding, and the reason π carries an informative prior at all.
- **The count intercept β₀ is deliberately non-identifiable** against the level of N.
  Expect its Rhat to be poor and say so; λ is unaffected.

### Results

1. Adult survival, and sex effect (males higher: OR 1.51, 0.99–2.32 — the interval includes 1; do not call it significant).
2. Tag effects on survival and on resighting.
3. Recruitment: delayed, split between ages 1 and 2, low overall.
4. Productivity.
5. Population trajectory and λ.
6. tLTRE: what sustains the population; what drives variation.

### Discussion

1. **Delayed recruitment** — natural-history contribution; consequences for any LEYE
   population model assuming age-1 breeding, including harvest assessments.
2. **Tag effects** — magnitude vs the geolocator literature; implications for
   deployment, permitting, and for inference in existing tracking studies.
3. **Adult survival in context** — relevance to harvest sustainability.
4. **What limits inference** — demographic stochasticity; recruitment/immigration
   confounding; relative-scale abundance. Frame the stochasticity result as a
   monitoring-design contribution, not an apology. Be explicit that the
   recruitment estimate is prior-driven while λ is not: λ holds at 0.929–0.931
   across all four prior settings — typical-cohort, cohort-averaged,
   1997-excluded and vague — while the recruitment share moves from 6.3% to 12.8%
   (see *Confounds to disclose*). Lead with that stability; it is what licenses the λ and tLTRE
   results even though the recruitment/immigration split is not identified.
5. **Conservation implications.**

---

## 6. Figures and tables

Most already exist at publication quality in `reports/leye_ipm_report.qmd`.

| # | Content | Status |
|---|---|---|
| F1 | Study area map | New |
| F2 | Life cycle with delayed recruitment and data sources | **Done** — `manuscript/figures/life_cycle.svg`. Colour-blind-safe (Okabe-Ito) with dashed recruitment arrows as a redundant cue. Cohort labels corrected to 122 in 1995–99 on 2026-09-30. |
| F2b | IPM as a directed acyclic graph | **Done** — `manuscript/figures/ipm_dag.svg`, in the Kéry & Schaub Fig 5.3 idiom. Its prior box reads "prior is the modern estimate" and names both cohorts; verified current 2026-10-22. |
| F3 | Population trajectory | Have |
| F4 | Adult survival by year | Have |
| F5 | Covariate forest plot (sex, tag; survival and resighting) | Have |
| F6 | tLTRE: level + variance decomposition | Have |
| F7 | Recruitment by era: modeled probability with individual cohorts, plus the posterior of the modern/1990s ratio | **Done** — `manuscript/figures/era_recruitment.png`, reproduced live in the report (§ *Has recruitment changed between eras?*) from the module posterior, so it cannot drift from the fit. Labels the 1997 cohort explicitly, since 7 of 11 returning carries most of the 1990s signal. |
| F7b | Same contrast, era estimates only | **Done** — `manuscript/figures/era_recruitment_simple.png`. No cohort scatter, so the y-axis closes in on the intervals; easier to read from the back of a room. Trade-off: the 1997 outlier is invisible, so pair it with the caveat in the text. Both variants come from the same script. |
| T1 | Priors and their sources | New |
| T2 | Posterior summaries of vital rates | Have (report tables) |
| T3 | Survival by sex × tag | Have |
| S1+ | Prior sensitivity, GOF, convergence | From Phase 1 |

---

## 7. Open questions

- Authorship and contributions — including the volunteer count programme and the
  1990s field crews.
- Data availability statement. The chick encounter history is now reproducible: cohort
  totals are real banding records held as fixed constants, and a rerun reproduces the
  committed file exactly (bar line endings). State that 106 of the 122 1990s rows are
  reconstructed never-resighted birds, and that brood IDs survive only for 1996.
- Whether to report the 1990s-era adult survival (`phi_pre` ≈ 0.51) as a historical
  comparison, or leave it as a nuisance parameter.
- How much modeling detail belongs in the main text vs supplement for this audience.

---

## References anchored so far

- Weiser et al. 2016, *Movement Ecology* — geolocator effects in 16 Arctic-breeding
  shorebirds; ~8%/yr reduction in apparent survival for leg-mounted tags.
- Brlík et al. 2020, *J. Animal Ecology* — meta-analysis; weak negative effect on
  apparent survival, scaling with relative load and harness attachment.
- Pakanen et al. 2020, *J. Avian Biology* — survival declines with time carrying the device.
- Koons et al. 2016, 2017 — transient LTRE.
- Dillenseger et al. 2026, *Oikos* — global shorebird survival; immature < adult,
  breeding-site immature estimates depressed by natal dispersal.
