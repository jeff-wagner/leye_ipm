# Manuscript plan

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
  modelling rests on BBS/CBC indices fitted with state-space logistic models — not
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
   by **two independent eras** — 139 chicks banded 1995–98 and 133 banded 2018–25 —
   with near-identical age-at-first-return patterns thirty years apart. Modern total
   probability of recruiting locally ≈ 5%.
2b. **Immigration, not local production, appears to sustain the population.** With the
   modern recruitment prior in place, local recruits account for only ~31% of annual
   gains (down from 56% under the 1990s-derived prior), and immigration contributes
   ~17% of λ against recruitment's ~7%. This is the headline the tLTRE now delivers,
   and it has direct conservation implications: a local population propped up by
   immigration is vulnerable to declines elsewhere in the flyway.
3. **GPS tags substantially reduce apparent survival** — a ~45% relative reduction,
   several times the published geolocator benchmark.
4. **Demographic stochasticity dominates variation in growth rate** in a population
   this small, which has direct consequences for monitoring design.
5. The population trend is negative but imprecise.

### Framing constraints — accept these early, do not overreach

- **λ = 0.967 (0.917–1.033). The interval includes 1.** This is *not* a demonstrated
  decline. Write "probably declining, most consistent with a slow decrease" and let
  the vital rates carry the paper.
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
   GPS. That is hard to explain as a modelling artifact.
2. **Not a detection artifact.** The model estimates tag effects on resighting
   separately. GPS birds show no detection penalty (OR 0.90, P(<1) = 0.59), so the
   survival effect is not birds becoming harder to see. Notably geolocator birds
   *do* show one (OR 0.52, P(<1) = 0.97) — worth reporting in its own right.
3. **Mechanistically expected direction.** Brlík et al. (2020) found effects scale
   with relative device load and with harness attachment; GPS tags increase both.

**Caveat that must be stated plainly:** this is *apparent* survival. Tag-induced
dispersal is not separable from mortality here.

---

## 2. Status

| Component | State |
|---|---|
| IPM structure, fitting, convergence | Done |
| tLTRE | Done |
| Quarto report with publication-quality figures | Done |
| Recruitment estimates | **Blocked — see Phase 0** |
| Goodness-of-fit / posterior predictive checks | **Missing** |
| Tag-effect robustness analyses | **Missing** |
| Novelty check (is this the first LEYE IPM?) | **Not done** |

---

## 3. Phase 0 — data resolution (blocking, do first)

**Downgraded but not resolved.** The IPM's recruitment prior now comes from the modern
era (133 chicks, 2018–2025), so the simulated 1990s cohorts no longer set it directly.
They still enter two ways: they size the 1990s half of the era comparison, and the 178
contemporaneous 1990s adults pin the shared detection that makes the modern estimate
identifiable. Five modern resights cannot determine detection alone.

So this is no longer blocking for the *headline* recruitment number, but it is still
blocking for any claim about **change between eras** — which is now a candidate result
(P(decline) = 0.91, ratio 0.10–1.61).

The 1990s cohort denominators for 1995/1997/1998 remain *simulated*
(`rpois(4, 0.85 * num_broods96)`, an explicitly arbitrary mean).

- [ ] Recover the real 1995/1997/1998 chick banding totals from the original records.
- [ ] Rebuild `data/LEYE_juv_EH_1995_2000.csv` from real denominators.
      `scripts/create_juv_eh_1995_2000.R` is now seeded, so the rebuild is reproducible.
      Delete the simulation branch once real totals are in.
- [ ] Confirm whether the 1990s chicks carried colour marks (the juvenile module
      currently shares detection between chicks and adults on that assumption; if they
      were metal-band-only, detection was lower and `psi` is underestimated).
- [ ] Refit `scripts/LEYE_juv90_module.R` → new `pi_rec` / `kappa` priors.
- [ ] Update the Beta priors in `scripts/LEYE_IPM.R` (currently `dbeta(2.268, 13.610)`
      and `dbeta(6.054, 6.309)`).
- [ ] Refit the IPM, re-run the tLTRE, re-render the report.
- [ ] **Re-check whether any conclusion changed.** Recruitment feeds λ and the tLTRE.

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
      **counts 0.323, productivity 0.512, mark-resight 0.473.**
- [x] **Prior sensitivity on recruitment.** Done, via prior hyperparameters passed
      as constants (so the sensitivity run is a config change, not a forked model).
      **Recruitment level is NOT prior-driven**: under a flat Beta(1,1), `pi_rec` =
      0.108 vs 0.106 informative. **Timing IS**: `kappa`'s interval widens from
      (0.23–0.75) to (0.02–0.97) with the mean unmoved. λ and the recruitment/
      immigration split are unchanged. Both fits converged (Rhat ≤ 1.004).
- [ ] **Novelty check.** Establish whether a field-based IPM exists for this species.
      Frame claims accordingly.

### Strengthening the tag result

- [x] **Time-since-tagging effect.** Done — `scripts/LEYE_tag_effects_CJS.R`.
      Duration slope −0.090 (−0.410, 0.225), P(<0) = 0.712: direction matches
      Pakanen but no detectable effect, and WAIC worsens (849.5 vs 845.9). The GPS
      effect survives its inclusion.
- [x] **Selection bias.** Done. **GPS birds were significantly heavier at banding:
      87.3 g vs 82.0 g untagged (+5.33, 95% CI 3.05–7.60, p < 0.001).** Geolocators
      were not (−0.20, p = 0.84). This works *against* the finding — heavier birds
      should survive better — so the tag effect is conservative.
- [x] **Standalone CJS check.** Done. M0 gives GPS OR 0.28 (0.15–0.53) vs the IPM's
      0.29 (0.15–0.54) — the tag effect is not an artifact of IPM linkage.
- [x] **Transience control.** Added because exposure differs sharply by tag class
      (untagged 25.2% first-intervals vs geo 9.0%, GPS 14.0%). Including it does not
      weaken the GPS effect (0.28 → 0.25).
- [ ] **Tag mass as % of body mass.** Use **87.3 g** as the GPS denominator, not the
      82.7 g population mean — GPS birds were a heavier subset. Still need tag specs.
- [ ] Attachment method (harness vs leg-mount) for each tag type — a documented
      effect modifier (Brlík et al. 2020).
- [ ] **Tag recovery dates**, if they exist. Geolocators must be recovered, so some
      birds stopped carrying tags mid-study; the duration test currently uses
      time-since-banding as an attenuated proxy. Recovery dates would sharpen it.

### Confounds to disclose in the manuscript

Found during the robustness work; all belong in the Methods or Discussion.

- **Tag type is largely confounded with banding era** — geolocators 2010–2018,
  GPS 2018–2020, untagged 1995–2025. Not fatal: 52 untagged birds were banded in
  the same 2018–2020 window as the GPS birds, so the contrast is not purely
  across eras. Year-varying baseline survival absorbs period effects.
- **Sex composition differs** — GPS 58.7% female vs untagged 44.8%. Since males
  survive better (OR 1.55), this biases GPS *downward*; it is controlled by the sex
  covariate, but the model must retain sex for that reason.
- **Transience runs the other way than expected** — the first-interval offset is
  *positive* (+0.659, 0.003–1.417), i.e. survival is higher just after banding.
  Mechanism unclear; worth a sentence rather than an explanation we cannot support.

### Optional but valuable

- [ ] Balanced-panel count sensitivity — `scripts/LEYE_trend_NB_coreSites_balancedPanel.R`
      already exists and addresses whether the trend is an artifact of growing site
      coverage (13 → ~32 wetlands).
- [ ] Power/monitoring framing: given demographic stochasticity dominates Var(λ),
      how long must this survey run to detect a decline of a given size? This turns
      a limitation into a contribution for a monitoring-oriented audience.

---

## 5. Phase 2 — structure

### Methods

1. **Study area and field methods** — Anchorage wetlands; volunteer survey design
   (3 weekly visits, 32 wetlands); banding and resighting; nest monitoring; tag
   deployment (types, masses, attachment, years).
2. **Data sources** — counts 2014–2025 (855 records); adult encounter histories
   2010–2026 (213 individuals); productivity 2018–2026 (77 nests); 1990s chick
   cohorts.
3. **IPM structure** — process model with delayed recruitment; the two-year lag and
   why; survival with sex and tag covariates; relative-scale count observation model.
4. **Recruitment sub-model** — the 1990s chick+adult joint analysis, marginalising
   over recruitment age. Explain *why*: a chick first seen at age 2 may have recruited
   at age 1 and been missed, so raw return rates confound timing with detection.
5. **tLTRE** — level and variance decomposition; explicit separation of demographic
   stochasticity.
6. **Implementation** — NIMBLE, chains/iterations/convergence, priors table.

### Results

1. Adult survival, and sex effect (males higher: OR 1.55, 1.01–2.39).
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
   monitoring-design contribution, not an apology.
5. **Conservation implications.**

---

## 6. Figures and tables

Most already exist at publication quality in `reports/leye_ipm_report.qmd`.

| # | Content | Status |
|---|---|---|
| F1 | Study area map | New |
| F2 | Life-cycle / model schematic with delayed recruitment | New — worth the effort, the structure is unusual |
| F3 | Population trajectory | Have |
| F4 | Adult survival by year | Have |
| F5 | Covariate forest plot (sex, tag; survival and resighting) | Have |
| F6 | tLTRE: level + variance decomposition | Have |
| T1 | Priors and their sources | New |
| T2 | Posterior summaries of vital rates | Have (report tables) |
| T3 | Survival by sex × tag | Have |
| S1+ | Prior sensitivity, GOF, convergence | From Phase 1 |

---

## 7. Open questions

- Authorship and contributions — including the volunteer count programme and the
  1990s field crews.
- Data availability statement. Note that the current chick encounter history predates
  the random seed and is not reproducible from the script; Phase 0 fixes this.
- Whether to report the 1990s-era adult survival (`phi_pre` ≈ 0.51) as a historical
  comparison, or leave it as a nuisance parameter.
- How much modelling detail belongs in the main text vs supplement for this audience.

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
