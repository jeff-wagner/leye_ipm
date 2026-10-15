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
   by **two independent eras** — 122 chicks banded 1995–99 and 133 banded 2018–25 —
   with near-identical age-at-first-return patterns thirty years apart. Modern total
   probability of recruiting locally ≈ 5%.
2b. **Immigration, not local production, appears to sustain the population.** With the
   modern recruitment prior in place, local recruits account for only ~32% of annual
   gains (down from 56% under the 1990s-derived prior), and immigration contributes
   ~17% of λ against recruitment's ~7%. This is the headline the tLTRE now delivers,
   and it has direct conservation implications: a local population propped up by
   immigration is vulnerable to declines elsewhere in the flyway.
   **State it with its range.** The ordering holds under both informative priors
   but the margin depends on which π the IPM receives: 16.7% vs 7.0% on the
   default, 14.3% vs 9.6% cohort-averaged — P(immigration > recruitment) 0.81 and
   0.67 respectively. Claim the direction, not the ratio. See *Confounds to
   disclose*.
3. **GPS tags substantially reduce apparent survival** — a ~43% relative reduction,
   several times the published geolocator benchmark.
4. **Demographic stochasticity dominates variation in growth rate** in a population
   this small, which has direct consequences for monitoring design.
5. The population trend is negative but imprecise.

### Framing constraints — accept these early, do not overreach

- **λ = 0.959 (0.914–1.015). The interval includes 1.** This is *not* a demonstrated
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
| Prior sensitivity | Done — and it **overturned** an earlier conclusion. Now three-way (typical-cohort, cohort-averaged, vague); λ is identical across all three, only the recruitment/immigration split moves |
| 2026 counts | **In** (added 2026-10-14). The series now runs 2014–2026, 935 records. λ 0.965 → 0.959 and P(λ<1) 0.88 → 0.94; the tLTRE split and every other conclusion held. Both sensitivity fits were re-run on the longer series |
| Model structure in the report | Done — life-cycle figure, DAG and the full statistical formulation are now a section of the Quarto report, matching *Model at a glance* below |
| Tag-effect robustness analyses | Done, but computed on the pre-correction encounter history — re-run |
| Recruitment estimates | Now from the **modern era**; the 1990s cohorts still enter via the era comparison and shared detection. Refitted 2026-09-29 on the real banding totals (122 chicks, 1995–1999) — no conclusion changed |
| Life-cycle figure and DAG | Done — both label the recruitment prior as the modern estimate, and name both cohorts (133 chicks 2018–25, 122 in 1995–99, with 178 1990s adults) |
| Novelty check (is this the first LEYE IPM?) | **Not done** |
| Study area map (F1), priors table (T1) | **Not done** |

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
- [ ] Confirm whether the 1990s chicks carried colour marks (the module shares detection
      between chicks and adults on that assumption; if they were metal-band-only,
      detection was lower and `psi` is underestimated). **Still open.**

**What moved.** Smaller 1990s denominators raise π₉₀ from 0.115 to 0.139, and because
detection is shared across eras that carries into the modern estimate too: π_M 0.049 →
0.054, κ 0.485 → 0.469. Downstream the changes are small — λ 0.963 → 0.965, the tLTRE
level split 76/17/7 → 75/17/7, P(immigration > recruitment) 0.83 → 0.81, GPS odds ratio
0.29 → 0.28. **Every conclusion in this document survives.** That robustness is itself
worth a sentence in the manuscript: the headline did not depend on the simulated
denominators, which is not something we could assert before.

**The one thing to disclose.** The 1997 cohort is 11 chicks with 7 returning — a 64%
local return rate against 11%, 7%, 0% and 0% in the other 1990s cohorts and 4% in the
modern era. 1997 was a genuinely low-effort season (21 adults marked, against
52/48/36/21 around it), so the small denominator is real, but one cohort now carries
most of the 1990s recruitment signal. Say so rather than leaving a reviewer to find it.

**Still blocking the era comparison**, though for a different reason than before: the
contrast is now P(decline) = 0.91 with a ratio spanning 0.08–1.59. That is suggestive,
not a result, and the width comes from 13 and 5 returning birds — not from the
denominators, which are now real.

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
      Now a **three-way** comparison — typical-cohort (default), cohort-averaged,
      and vague — reported in the Quarto report's prior-sensitivity section with
      π, κ, ω, λ and the contribution shares for each. Both companion fits are
      built by `scripts/LEYE_IPM_sensitivity.R` (no argument runs both;
      `cohortmean` or `vague` runs one), which takes the data, structure and
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
      0.28 (0.15–0.53) — the tag effect is not an artifact of IPM linkage.
      *Caveat:* the standalone fit is from 2026-08-06, before the encounter-history
      correction, so it is the one side of this comparison still to be re-run
      (see the Status table). The IPM side is current as of 2026-10-14.
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
- **Recruitment is prior-driven, and immigration absorbs the difference**
  (re-checked 2026-10-02 across three prior settings). Under the default prior
  the posterior barely leaves it: π moves 0.0536 → 0.0513, **−0.0023**, and the
  interval is essentially unchanged (prior 0.006–0.149, posterior 0.006–0.135).
  The counts carry almost no independent recruitment signal, and none at all
  about the age split — under vague priors κ returns 0.491 (0.024–0.975), which
  *is* Beta(1,1).

  | prior on π | π | ω | λ | rec. share | imm. share | P(imm > rec) |
  |---|---|---|---|---|---|---|
  | typical cohort (default) | 0.051 | 0.157 | 0.961 | 7.0% | 16.7% | 0.81 |
  | cohort-averaged | 0.071 | 0.135 | 0.961 | 9.6% | 14.3% | 0.67 |
  | vague | 0.103 | 0.114 | 0.963 | 13.5% | 12.1% | 0.46 |

  The three runs trace the recruitment–immigration ridge: as π rises, ω falls to
  compensate and **λ does not move** — identical to three decimals. The paper's
  growth-rate results are immune to this choice; only the *attribution* is not.
  Disclose it rather than let a reviewer find it: it is the empirical signature
  of the additive ψ/ω confounding set out in *Model at a glance*, and it is
  precisely why π carries an informative prior. The honest framing is that the
  chick data identify recruitment and the counts identify λ; neither alone
  separates recruitment from immigration.

- **Which π the IPM should receive is a genuine modelling choice** (added
  2026-10-02). `LEYE_recruitment_two_era.R` fits a cohort random effect on the
  logit scale and reports π_M = ilogit(μ_π + b_era) — the *typical* cohort, with
  its random effect at zero. But the IPM uses π in an expected-count expression,
  R = N·(f/2)·ψ, and has no cohort random effect of its own, so arguably it wants
  E[π_c], the mean **across** cohorts. With σ_π ≈ 1.39 the two differ by about a
  third (0.054 vs 0.072) through Jensen's inequality. Both are now fitted and
  reported side by side; neither is obviously wrong. Hold the cohort-averaged
  version loosely, because σ_π rests on 13 cohorts and 18 returning birds and the
  1997 cohort alone sits at π_c = 0.58 against 0.03–0.18 elsewhere — and since
  σ_π is shared across eras, that one cohort inflates the gap for the modern
  cohorts too. **Report the immigration ordering at the conservative end:**
  immigration exceeds local recruitment with P ≈ 0.67–0.81 depending on the
  definition, a margin of roughly 1.5:1 to 2.4:1 rather than a settled number.

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
2. **Data sources** — counts 2014–2026 (935 records); adult encounter histories
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
   monitoring-design contribution, not an apology. Be explicit that the
   recruitment estimate is prior-driven while λ is not: λ holds at 0.959–0.963
   across all three prior settings — typical-cohort, cohort-averaged and vague —
   while the recruitment share moves from 7% to 13.5% (see *Confounds to
   disclose*). Lead with that stability; it is what licenses the λ and tLTRE
   results even though the recruitment/immigration split is not identified.
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
