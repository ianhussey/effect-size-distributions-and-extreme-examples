# Project context: Beta-bound SMD plausibility

## What this project is about

Ian (a quantitative psychology researcher) is developing a principled method for detecting implausibly large standardised mean difference (SMD) effect sizes on bounded psychological scales — the kind of scales that are typically sum scores of Likert items, with a logical floor and ceiling.

The core observation that motivates the project: many psychological scales have a logical range that is roughly 4–9 times their SD. A naive linear bound says that if the mean sits at, say, 75% of the maximum and the range is 4 SD, you cannot observe an upward shift of more than 1 SMD or a downward shift of more than 3 SMD before bumping into the scale's bounds. But this naive bound implicitly assumes the SD is a fixed constant. In reality, when the mean of a bounded variable approaches a bound, the SD must compress — there is less room for the distribution to spread. Ignoring this compression overstates how large an SMD can plausibly be.

The project's premise is that fitting a Beta distribution to the observed (M, SD) on the rescaled (0, 1) interval, and then propagating the implied precision parameter under a constant-precision assumption, gives a tighter and more honest ceiling on plausible SMDs than the naive linear bound. This puts the work in the same family as GRIM and SPRITE-style plausibility checks but targets SMDs rather than (M, SD) pairs directly.

## The math, briefly

Rescale the scale to (0, 1) so that mu = (M − min) / (max − min) and sigma_tilde = SD / (max − min). For a Beta distribution with mean mu and precision phi, variance is mu(1 − mu) / (1 + phi). Method of moments inverts this: phi = mu(1 − mu) / sigma_tilde^2 − 1, alpha = mu * phi, beta = (1 − mu) * phi.

Assuming phi is constant across conditions (the load-bearing assumption), the SD predicted at any other mean mu' is sqrt(mu'(1 − mu') / (1 + phi)) * range. This curve is what lets you ask "how many SDs of change are possible before hitting the bound?" in a way that respects the geometry of bounded scales.

Three things to keep distinct when computing SMD:

1. **Naive (control SD)**: SMD = (M' − M) / SD_obs. This is identical to the linear bound; Beta regression doesn't change it.
2. **Pooled SD (Beta-predicted)**: SMD = (M' − M) / sqrt((SD_obs^2 + SD_post^2) / 2). Bounded; in the limit mu' → bound, asymptotes to sqrt(2) × the naive bound.
3. **Post-only SD (Beta-predicted)**: SMD = (M' − M) / SD_post. Diverges to infinity near the bounds because SD_post → 0. Generally a poor convention for plausibility checks.

The Bhatia–Davis inequality gives a parameter-free sanity floor: for any distribution on [min, max] with mean M, SD ≤ sqrt((M − min)(max − M)). A reported SD that exceeds this is mathematically impossible regardless of distributional assumptions. Beta regression provides a tighter, parametric refinement of the same geometry.

## What's in this folder

- `beta_smd_explorer.html` — standalone single-file HTML/JavaScript app. Open in any browser. Lets you change scale bounds, M, and SD via inputs and sliders, and shows the implied Beta PDF, the SD-vs-mean curve with the Bhatia–Davis envelope overlaid, the SMD-vs-target-mean curves under the three denominator conventions, and the limiting SMD ceilings. The original prototype.
- `beta_smd_bounds.R` — clean R script with the core functions (`bhatia_davis_max_sd`, `beta_fit_mom`, `predict_sd_at_mean`, `smd_bounds`, `smd_curves`) plus a base-R `plot_beta_smd` demo. Designed to be sourced and applied row-wise to meta-analytic datasets to flag suspect SMDs.
- `beta_smd_shiny/app.R` — Shiny app version of the HTML explorer (uses ggplot2). Run with `shiny::runApp("beta_smd_shiny")` or click "Run App" in RStudio.
- `beta_smd_app.R` — earlier copy of the Shiny app code at the top level of the folder; superseded by `beta_smd_shiny/app.R` (RStudio's Run App button needs the standard `app.R`-in-folder layout).

## Open questions and known caveats

The constant-phi assumption is the big one. Real interventions can shift variance independently of the mean, in which case the Beta-implied SD curve is misleading. The right framing is that the curve produces a *reference distribution* of SMDs consistent with a constant-precision Beta data-generating process, not a hard cap.

Likert sum scores are discrete and often multimodal in ways Beta won't capture. Zero-and-one-inflated Beta extensions exist for scales with substantial floor/ceiling mass and may be worth pursuing.

The next concrete empirical step is to simulate or fit this against real Likert scales to see how tight the Beta bound actually is in practice — i.e. whether it earns its keep over the back-of-envelope linear bound for realistic effect sizes encountered in published meta-analyses.

## How Ian likes to be engaged on this

He's methodologically literate and wants the math at full depth, with assumptions surfaced and caveats stated. Don't oversell — be explicit about what an approach buys you and what it doesn't (especially the parametric vs. non-parametric distinction between Bhatia–Davis and the Beta bound).
