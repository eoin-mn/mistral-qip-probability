# Statistical model — EmergApps end-of-project assessment (v11)

A one-page technical summary of the model implemented in `end-of-project-assessment.html`. Intended for a statistical reviewer. The tool is an **end-of-project Bayesian assessment** for Emergency Department QI projects: it answers two graded questions — *did the process genuinely improve?* and *is the true gain as big as a pre-specified target?* — rather than returning a p-value or a binary verdict.

## 1. Data and measurement model

The user enters raw charted points grouped into phases: a baseline phase (B points) and one to three post-change PDSA cycles. The judged "end" set is either **all points since the first change (pooled)** or **only the final cycle's points** (user's choice, declared before assessment). Every charted point is the **median of n patient-level observations**, n declared by the user.

- **Continuous data (e.g. times).** Patient-level values are lognormal: in a phase with true median level m, a patient's value is X = m·exp(σZ), Z ~ N(0,1). A charted point is the median of n such values. The simulation samples the median of n normals **exactly from its order-statistic distribution** (odd n: Φ⁻¹(U), U ~ Beta((n+1)/2, (n+1)/2); even n: the average of the two central order statistics, sampled from their joint distribution), so small-n median uncertainty is modelled exactly, not by a normal approximation.
- **Proportions.** A charted point is K/n with K ~ Binomial(n, p) (normal approximation only when n·p·(1−p) ≥ 25). Entered percentages are scaled to (0,1); fractional implied counts (e.g. 75% from 6 patients = 4.5 patients) are **kept fractional**, never rounded.
- **σ** is either chosen from three presets (log-scale SD 0.195 / 0.514 / 0.93) or estimated from the baseline by successive log-differences,

  σ̂² = [ Σᵢ (log xᵢ₊₁ − log xᵢ)² / (2(B−1)) ] / v(n),  with v(n) = π/(2n),

  where π/(2n) is the asymptotic variance of a median of n normals; dividing it out converts point-median scatter back to patient-level spread. The factor is slightly conservative at small n.

## 2. Effect measure and hypotheses

Let δ denote the true improvement in the beneficial direction:

- continuous: **relative** change in the true median level (δ = (m_post − m_base)/m_base, signed for benefit);
- proportions: **absolute** change in the true rate (percentage points).

The observed shift t_obs is the same functional applied to the charted medians. Two hypotheses are compared with equal prior weight:

- H₀ (spike): δ = 0, prior π₀ = ½;
- H₁ (slab): δ ~ Uniform(0, d_max], prior π₁ = ½, where d_max = max(5·t\*, 1.5·|t_obs|) and t\* is the user's pre-specified target.

Numerically the slab is discretised into 40 **geometric** cells on (t\*/10, d_max], cell c carrying weight w_c = hi − lo and representative δ_c = √(lo·hi).

## 3. Likelihood: the tail event

For any candidate δ, the engine simulates R paired pseudo-datasets — a baseline median and an end median under true effect δ — using the exact small-n samplers above (baseline draws are reused across δ: common random numbers). Let T be the resulting shift statistic. The evidence is coarsened to the **tail event** E = {T ≥ t_obs − ε}, with likelihood

  ℓ(δ) = log P(E | δ) ≈ log [ (#{r : T_r ≥ t_obs − ε} + ½) / (R + 1) ],

the add-half smoothed tail frequency. Under H₀ the same quantity ℓ(0) is the chart-luck (false-alarm) probability.

**Why coarsen?** With heavily quantized data (e.g. ten consecutive identical medians) a pointwise likelihood is degenerate or unstable, while the tail event "a chart at least this good" remains well defined and monotone in δ. The cost is loss of fine information; the benefit is robustness. This choice is what fixed early under-confidence on quantized datasets.

## 4. Posterior quantities

**Chance the improvement is real.** With log m₁ = logSumExp_c[ℓ(δ_c) + log w_c] − log d_max,

  P(δ > 0 | data) = 1 / (1 + exp(−(log m₁ − log m₀))),   m₀ = exp ℓ(0).

**Chance the true gain reaches the target.** P(δ ≥ t\* | data) = P(δ > 0 | data) × P(δ ≥ t\* | δ > 0, data).

- *Proportions — exact branch.* With K = Σ x_j·n and N = n·P, a flat prior gives the true end rate p₁ | data ~ **Beta(1 + K, 1 + N − K)**. The target is anchored on the **measured** baseline: threshold t = p̂₀ + t\* (or p̂₀ − t\* if lower is better), and P(δ ≥ t\*) = P(p₁ ≥ t), evaluated as 1 − I_t(a, b) via the regularised incomplete beta function (Lanczos log Γ plus Lentz continued fraction, valid for fractional a, b). Effect estimates are posterior quantiles of p₁ − p̂₀ (10/50/90th). A useful calibration property: the posterior median of p₁ equals the entered end level.
- *Continuous — tail-curve branch.* The points (δ, P(E|δ)) over the grid plus δ = 0 and δ = t\* are made monotone non-decreasing; the effect estimate is δ̂ = inf{δ : P(E|δ) ≥ ½} — "which true effect would make a chart like yours a typical, rather than lucky, outcome" — with an 80% interval from the 0.1 and 0.9 crossings; P(δ ≥ t\*) = 1 − P(E | t\*).

**Luck envelope.** The 2.5th–97.5th percentiles of the end median under H₀, displayed as "if your change had done nothing, your end median would plausibly have landed between …".

**Wrong-way data.** If t_obs ≤ 0 the run is mirrored (direction reversed, |t_obs|) and the display reports 1 − P(real); the target question is marked not applicable.

R = 2,500 headline draws (800 per cycle-table row); figures wobble by roughly ±1 percentage point. Computation is chunked to keep the page responsive.

## 5. Data-quality checks (advisory, not model inputs)

Baseline steadiness: a Theil-slope trend test (|slope|·√(S_xx) against a successive-difference noise estimate σ̂_d = √(ΣΔ²/(2(B−1))), flagged at z > 2.5); a first-half vs second-half median comparison; a repeated-values check (>50% identical points flagged, as repeats make medians jump). Sample-size advisories: ≥8 baseline points, ≥4 end points, ≥5 patients per point.

## 6. Assumptions

1. Charted points are **independent**; no autocorrelation, trend, or seasonality within any phase.
2. Every point rests on the same declared n; patient-level values are lognormal (continuous) or Bernoulli with common p (proportions).
3. The change takes effect immediately and is constant within each cycle; σ is common to all phases.
4. Nothing else changed contemporaneously — surges, staffing, case mix and co-interventions are unmodelled (no covariate adjustment).
5. Prior framing: 50:50 spike-vs-slab, slab uniform to d_max. Results are conditional on this choice.
6. In the proportion target branch, the baseline is anchored at its **measured** level p̂₀ (treated as known); the headline branch does model baseline uncertainty.
7. The target t\* is honest and pre-specified; the tool cannot verify this.

## 7. Strengths

- **Graded Bayesian answers** aligned to the decision ("real?", "target reached?") rather than a threshold ritual; the target is declared before looking at the result.
- **Exact small-n geometry**: medians-of-n are simulated from their true order-statistic distributions, so thin data are penalised naturally and correctly.
- **Robust to quantized/repeated data** via the tail-event likelihood, where pointwise likelihoods collapse.
- **Winner's-curse control**: the effect estimate is the δ that makes the observed chart *typical*, not the observed shift itself; the spike-at-zero prior guards against declaring luck.
- **Exact beta-binomial branch** for proportions, with fractional counts preserved and closed-form target/effect figures.
- **Empirically calibrated** (200-run simulation checks): when the true effect equals the target, P(reached) ≈ 38–54%; above it, 98–100%; below it, ≈0–1%; false-alarm rate at the 95%"real" line ≈ 1.3–4.7% across scenarios.
- Self-contained single file; all assumptions disclosed in-page; no data leaves the browser.

## 8. Weaknesses and caveats

- **Independence is the big one.** ED series drift and cluster (autocorrelation, special causes); the model sees pseudo-replication, so probabilities are **flattering** when points move together. This is disclosed in-page but not corrected.
- **No adjustment for concurrent events** (surges, staffing, seasonality, case mix). A shift that coincides with such an event will be attributed to the change.
- **Prior sensitivity.** With weak data, p(real) depends materially on the uniform-slab ceiling d_max; the 50:50 starting point is a fixed, debatable choice.
- **Monte Carlo wobble** of about ±1 pp; small-cell discretisation of the slab; the continuous target branch (tail-curve inversion) is coarser than the exact proportion branch.
- **Multiple comparisons**: the cycle-by-cycle table runs several unadjusted tests; a single row can look convincing by chance. Users are told to weight the headline.
- **Baseline anchoring in the target branch** ignores baseline uncertainty for proportions; a mis-measured baseline shifts the threshold.
- σ̂ from a short baseline is itself uncertain (floor of 0.05 imposed); preset σ values are assumptions, not measurements.
- **Not an SPC monitor**: it is a retrospective, whole-project assessment. It does not detect drift, and run-chart/SPC rules remain the right tool for prospective monitoring.
