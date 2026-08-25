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

- **λ = 0.963 (0.914–1.025). The interval includes 1.** This is *not* a demonstrated
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

> **Numbers in this document are a snapshot, and the fit moves.** Everything below was
> checked against `LEYE_IPM_samples.rds` on 2026-08-14. λ, the goodness-of-fit
> p-values and the sex effect all drifted slightly when the modern recruitment prior
> went in — small shifts, but the sex effect's lower bound moved onto exactly 1.00,
> which changes how it can be described. Re-verify every quoted value against the
> current posterior before drafting, and again before submission.

| Component | State |
|---|---|
| IPM structure, fitting, convergence | Done |
| tLTRE | Done |
| Quarto report with publication-quality figures | Done |
| Goodness-of-fit / posterior predictive checks | Done — and it caught a real misspecification |
| Prior sensitivity | Done — and it **overturned** an earlier conclusion |
| Tag-effect robustness analyses | Done, but computed on the pre-correction encounter history — re-run |
| Recruitment estimates | Now from the **modern era**; the 1990s cohorts still enter via the era comparison and shared detection |
| Life-cycle figure and DAG | Done — both carry a stale "1990s prior" label |
| Novelty check (is this the first LEYE IPM?) | **Not done** |
| Study area map (F1), priors table (T1) | **Not done** |

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
- [ ] Refit `scripts/LEYE_recruitment_two_era.R` → new `pi_M` / `kappa` priors.
      (Superseded `LEYE_juv90_module.R`, which fitted the 1990s alone; that script
      is retained only for its validated likelihood functions.)
- [ ] Update the Beta priors in `scripts/LEYE_IPM.R` (currently `dbeta(2.4200, 47.2852)`
      and `dbeta(7.7764, 8.2649)`, from the modern era).
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
      **counts 0.313, productivity 0.498, mark-resight 0.475.**
- [x] **Prior sensitivity on recruitment.** Done, via prior hyperparameters passed
      as constants (so the sensitivity run is a config change, not a forked model).
      **This conclusion REVERSED when the modern recruitment prior replaced the
      1990s one, and the earlier version of this entry was wrong.**

      Under the 1990s prior, informative and vague agreed (0.106 vs 0.108) and the
      recruitment level looked safely data-driven. Under the modern prior they
      differ two-fold — `pi_rec` = **0.050 informative vs 0.105 vague** — because
      the counts and the modern chick resights disagree, and the tighter modern
      prior wins. The old agreement was a coincidence: the 1990s prior happened to
      sit near what the counts wanted.

      **This belongs in the Discussion, not buried in a supplement.** Three readings
      cannot be separated with current data: (a) the counts cannot distinguish
      recruitment from immigration, so count-inferred recruitment absorbs
      unmodelled inflow; (b) assuming modern chick detection matches the 1990s may
      be too generous, biasing recruitment down; (c) recruitment genuinely is low
      and immigration genuinely high.

      **Timing** is still prior-driven: `kappa`'s interval widens (0.25–0.72) to
      (0.02–0.97) with its centre unmoved. Both fits converged (Rhat ≤ 1.004).
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
- [ ] **Re-run the tag analyses against the corrected encounter history.** Every tag
      number quoted below (GPS OR 0.28–0.29, the transience and duration models, the
      morphometric comparisons) was computed on the 213-bird history, before band
      129232054 was removed and 1422220945 corrected. The change is one bird out of
      213 and neither is a tagged bird, so the estimates should barely move — but
      "should barely move" is an assumption, and the published numbers need to come
      from the data actually in the repository. Re-run `LEYE_tag_effects_CJS.R` and
      the selection-bias check before drafting.
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
  survive better (OR 1.53), this biases GPS *downward*; it is controlled by the sex
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
   2010–2026 (**212** individuals); productivity 2018–2026 (77 nests); chick
   cohorts from **two eras** — 139 banded 1995–98 and 133 banded 2018–25.

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

Transcribed from `scripts/LEYE_IPM.R` on 2026-08-14. Indices: *t* = year (1…T, 2014–2025),
*i* = count record, *j* = marked individual, *m* = nest, *w* = survey week, *k* = wetland.
**If the model changes, this block must be re-transcribed — it is documentation, not code.**

```
──── PROCESS ─────────────────────────────────────────────────────────────────
  N₁, N₂        ~  Normal(0, 50) truncated at 0          two free initial years,
                                                          because recruitment lags 2 yr
  for t = 3…T:
    Rₜ          ~  Poisson( N_{t-1}·(f_{t-1}/2)·ψ₁  +  N_{t-2}·(f_{t-2}/2)·ψ₂ )
    Sₜ          ~  Poisson( N_{t-1}·φ_{t-1} )
    Iₜ          ~  Poisson( N_{t-1}·ω )
    Nₜ          =  Sₜ + Rₜ + Iₜ                          recruits enter the breeding
                                                          class directly; no juvenile stage
    λₜ          =  N_{t+1} / Nₜ

──── VITAL RATES ─────────────────────────────────────────────────────────────
  logit φₜ      =  μ_φ + ε_φ,ₜ            ε_φ,ₜ ~ Normal(0, σ_φ)   adult survival
  logit sₜ      =  μ_s + ε_s,ₜ            ε_s,ₜ ~ Normal(0, σ_s)   nest success
  fₜ            =  sₜ · C · q                         C = 4 (clutch), q = per-egg fledging
  ψ₁            =  π · κ                              recruit at age 1
  ψ₂            =  π · (1 − κ)                        recruit at age 2
  ω             =  immigration rate (constant — see the note at its prior)

──── OBSERVATION 1 · COUNTS  (relative index) ────────────────────────────────
  log μᵢ        =  β₀ + log N_{t(i)} + week_{w(i)} + b_e·E1ᵢ + b_d·E2ᵢ
                     + u_{k(i)} + δᵢ
  countᵢ        ~  Poisson(μᵢ)
      u_k ~ Normal(0, σ_site)      wetland random effect
      δᵢ  ~ Normal(0, σ_od)        overdispersion
      β₀ is non-identifiable against the level of N — read shape, not magnitude

──── OBSERVATION 2 · MARK–RESIGHT  (CJS, latent alive-state z) ───────────────
  logit φ_{j,t} =  ℓφ[ lookup(t) ] + βˢᵉˣ_{sex(j)} + βᵗᵃᵍ_{tag(j)}
  logit p_j     =  p₀ + βᵖ_{tag(j)}
  z_{j,t}       ~  Bernoulli( φ_{j,t-1} · z_{j,t-1} ),   z_{j,f(j)} = 1
  y_{j,t}       ~  Bernoulli( p_j · z_{j,t} )
      ℓφ[1…T-1] = logit φₜ ;  ℓφ[T] = logit φ_pre ;  ℓφ[T+1] = logit φ_post
      the two overflow slots absorb intervals outside the count series

──── OBSERVATION 3 · PRODUCTIVITY  (two-part) ────────────────────────────────
  succ_m        ~  Bernoulli( s_{t(m)} )               did the nest fledge anything?
  chicks_m      ~  Binomial( C, q · succ_m )           how many of the clutch
      a single Poisson was REJECTED here: p = 0.000, because the data are
      bimodal (0, 2, 3, 4 chicks; none with exactly 1)

──── PRIORS ──────────────────────────────────────────────────────────────────
  μ_φ, μ_s, p₀, βˢᵉˣ, βᵗᵃᵍ, βᵖ  ~ Normal(0, 1.5)
  β₀, b_e, b_d, week_w          ~ Normal(0, 5)
  σ_φ ~ Exp(2)   σ_s, σ_site, σ_od ~ Exp(1)   ω ~ Exp(1)
  φ_pre, φ_post, q              ~ Uniform(0, 1)
  π   ~ Beta(2.4200, 47.2852)   ← modern-era recruitment, external module
  κ   ~ Beta(7.7764,  8.2649)   ← age-1 fraction, external module
      these two are the ONLY informative priors, and both are passed as
      constants so the vague-prior sensitivity run is a config change
```

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

1. Adult survival, and sex effect (males higher: OR 1.53, 1.00–2.34 — lower bound sits on 1; do not call it significant).
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
| F2 | Life cycle with delayed recruitment and data sources | **Done** — `manuscript/figures/life_cycle.svg`. Colour-blind-safe (Okabe-Ito) with dashed recruitment arrows as a redundant cue. **Legend still says the recruitment prior is from the 1990s cohorts; it is now modern-era — needs updating.** |
| F2b | IPM as a directed acyclic graph | **Done** — `manuscript/figures/ipm_dag.svg`, in the Kéry & Schaub Fig 5.3 idiom. Same stale "1990s chick cohorts / informative prior" box as F2. |
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
