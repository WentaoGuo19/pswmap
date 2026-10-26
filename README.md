# pswmap

**pswmap** borrows external control data, such as historical trials or real-world data, for the control arm of a
randomized trial. It does so through a **propensity-score-weighted meta-analytic predictive (PSW-MAP) prior**, and it
handles baseline covariates that are **missing in the external data**.

It implements the method of Guo and Zhu, *Propensity Score Integrated MAP Prior with Historical Missingness*. The
simulation study of the paper can be reproduced with the companion repository
[pswmap-paper](https://github.com/WentaoGuo19/pswmap-paper).

## What the package does

You have a randomized trial with treated and control subjects, and one or more external sources of control subjects.
The external subjects differ from the trial population in their covariates, and some of their covariates may be
missing. The trial is complete. `pswmap` estimates the treatment effect in four steps:

| Step | Function | What it does |
|---|---|---|
| 1 | `psw_weights()` | For each external source, fits a propensity score to the trial stacked with that source only, and forms the odds weights $w = e/(1-e)$. These reweight the source to the trial population. |
| 2 | `psw_summaries()` | Summarizes each source by its weighted mean and design-based variance. Under multiple imputation, the summaries are combined across imputed datasets by Rubin's rules. |
| 3 | `psw_map_prior()` | Combines the source summaries in a normal-normal hierarchical model and returns the robust MAP prior for the trial's control mean. Sources with poor overlap have a larger variance and are discounted automatically. |
| 4 | `psw_posterior()` | Updates the priors with the trial data and returns the posterior of the treatment effect. |

`psw_reference_scale()` computes the reference scale $\sigma_{ref}$ used by steps 3 and 4.

`pswmap()` runs all four steps in one call and is the function most users need. With `borrow = FALSE` it analyses
the trial alone, which gives the reference analysis without borrowing.

## How the method works

Let $Z = 1$ mark trial subjects and $Z = 0$ external subjects, and let $H_1, \dots, H_K$ be the external sources.

**1. Weights.** External subjects are not randomized, so their covariates $X$ follow a different distribution
$f_k(x)$ than the trial's $f_1(x)$. For each source, the propensity score $e_k(x) = P(Z = 1 \mid X = x)$ is fitted
to the trial stacked with that source only. Its odds are proportional to the density ratio $f_1(x)/f_k(x)$, so the
odds weight

$$w_i = \frac{e_k(x_i)}{1 - e_k(x_i)}$$

reweights the source to the trial population. The score is capped at 0.95, so no weight exceeds 19.

**2. Source summaries.** Each source is summarized by its weighted (Hajek) mean and design-based variance,

$$\bar y_k^w = \frac{\sum_i w_i y_i}{\sum_i w_i}, \qquad V_k = \frac{\sum_i w_i^2 (y_i - \bar y_k^w)^2}{(\sum_i w_i)^2}.$$

**3. MAP prior.** The summaries enter a normal-normal hierarchical model (Neuenschwander et al., 2010),

$$\bar y_k^w \mid \theta_k \sim N(\theta_k, V_k), \qquad \theta_0, \theta_1, \dots, \theta_K \mid \mu, \tau \sim N(\mu, \tau^2), \qquad \tau \sim \text{half-}N\big(0, (\sigma_{ref}/c)^2\big),$$

where $\theta_0$ is the trial's control mean. Its predictive distribution is the PSW-MAP prior. Sources that agree
give a small $\tau$ and much borrowing, and a source with poor overlap (large $V_k$) is discounted. A weakly
informative component of weight 0.2 is added (Schmidli et al., 2014). It takes over when the trial data conflict
with the external sources.

**4. Posterior.** The prior is updated with the trial's control arm, the treated arm gets a weakly informative
prior, and the posterior of the treatment effect $\beta_1 = \theta_T - \theta_0$ is reported, with the decision
rule $P(\beta_1 > 0 \mid \text{data}) > 0.975$.

The method assumes that, given the measured covariates, external and trial controls have the same outcome
distribution, apart from the between-source differences that $\tau$ and the robust component absorb. It also
assumes that every source has subjects resembling the trial subjects (overlap). The help pages (`?pswmap`,
`?psw_weights`, `?psw_summaries`, `?psw_map_prior`, `?psw_posterior`) give the theory of each step with formulas
and references.

## Installation

```r
# install.packages("remotes")
remotes::install_github("WentaoGuo19/pswmap")
```

The default analysis needs only RBesT and mice, which are installed with the package. The flexible propensity score
learners need extra packages: bartMachine for `"bart"` (which also needs Java), ranger for `"rf"` and xgboost for
`"xgb"`.

## Quick start

```r
library(pswmap)
data(pswmap_example)
trial    = pswmap_example$trial      # outcome Y, treatment trt (1/0), covariates X1-X10
external = pswmap_example$external   # list of three sources (H1, H2, H3) with Y and X1-X10, some covariates missing

fit = pswmap(trial, external, tau_c = 16, seed = 1)
fit
```

```
PSW-MAP analysis: 3 external source(s), learner logit | missing covariates: mi (5 imputations)
Heterogeneity constant c = 16 | reference scale = 1.3 | prior ESS = 142

External sources:
 source   n n_incomplete mean    se kish_ess relweight_incomplete
     H1 200           47 2.03 0.102    175.5                0.914
     H2 200           48 2.01 0.113     81.7                1.105
     H3 200           33 1.88 0.106    156.7                0.910

Treatment effect (treated - control):
  posterior mean 0.769 (SD 0.147) 95% CrI [0.491, 1.07]
  P(effect > 0) = 1 -> effective at threshold 0.975
```

The output contains:

- **External sources**: one row per source. It shows the weighted mean of the outcome and its standard error, the
  Kish effective sample size of the weights, and the relative weight of incomplete subjects. The relative weight is 0
  if incomplete subjects are dropped.
- **Prior ESS**: how many trial control subjects the robust prior is worth.
- **Treatment effect**: the posterior of the treated-minus-control difference in means, and the probability that it
  is positive.

`summary(fit)` also prints the prior itself.

## Choosing the options

**Data format.**
- `trial` is a data frame with the outcome, a 0/1 treatment column and the covariates.
- `external` is either a named list of data frames, one per source, or one data frame with a `source` column.
- Covariates must be numeric: code binary covariates 0/1, and expand factors into indicator columns.
- By default every column other than the outcome and treatment is used as a covariate; use `covariates` to choose.
- The outcome must be continuous.

**Missing external covariates** (`missing`). Covariates are *missing completely at random* (MCAR) when
missingness is independent of the data, *missing at random* (MAR) when it depends only on observed values (other
covariates or the outcome), and *missing not at random* (MNAR) otherwise (Rubin, 1976).

| Option | What it does | Valid when |
|---|---|---|
| `"mi"` (default) | multiple imputation within each source, with the outcome in the imputation model, combined by Rubin's rules | MAR, including missingness that depends on the outcome |
| `"cc"` | complete-case weighting: incomplete external subjects get weight 0 | missingness independent of the outcome given the covariates (borrows less) |
| `"native"` | the learner's own missing-data handling (missingness indicators, BARTm, random-forest or XGBoost split directions) | as complete-case: because the trial is complete, a missing value marks a subject as external, so these methods move toward complete-case weighting |

Multiple imputation is the recommended default. When missingness depends on the outcome, complete-case weighting and
native handling bias every source in the same direction. The robust MAP prior cannot detect that bias, because the
sources still agree with one another. Control the imputation with `mi_m` (number of imputed datasets, default 5),
`mi_engine` (`"default"`, `"cart"` or `"rf"`) and `mi_args` (passed to `mice::mice()`).

**Propensity score learner** (`learner`).
- `"logit"` (default): main-effects logistic regression.
- `"bart"`, `"rf"`, `"xgb"`: Bayesian additive regression trees, random forests and XGBoost.

In the simulation study, logistic regression gave the lowest mean squared error, even when the true propensity
score was nonlinear. Tune the flexible learners with `learner_args`. The propensity score is capped at
`ps_cap = 0.95` before the weight is formed, so no external subject gets a weight above 19.

**Heterogeneity constant** (`tau_c`, required). The between-source standard deviation $\tau$ has a half-normal prior
with scale $\sigma_{ref}/c$. Use `tau_c` = 32, 16, 8 or 4 for small, moderate, substantial or large anticipated
heterogeneity. A smaller value borrows less and gives more protection against prior-data conflict. To see how much
the conclusion depends on this choice, rerun the analysis with the same `seed`. The propensity score weights are then
identical, and only the prior on $\tau$ changes:

```r
fits = lapply(c(32, 16, 8, 4), function(c) pswmap(trial, external, tau_c = c, seed = 1))
sapply(fits, function(f) f$effect["mean"])
```

**Robust prior and posterior.**
- `robust_weight` (default 0.2), `robust_mean` and `robust_sd` set the weakly informative component.
- `treated_prior` sets the prior of the treated arm.
- `threshold` (default 0.975) sets the decision rule $P(\text{effect} > 0) > $ `threshold`.
- `map_control()` holds the MCMC and mixture settings of the MAP prior.

**The trial alone.** `pswmap(trial, borrow = FALSE)` analyses the trial without external data. It is the natural
benchmark for the gain from borrowing.

## Step by step (for advanced users)

`pswmap()` is all that most analyses need. If you want to compute and inspect each step yourself, the four steps are
also available as separate functions:

```r
set.seed(1)
w     = psw_weights(trial, external)                       # 1. propensity scores and odds weights
s     = psw_summaries(w)                                   # 2. weighted means, SEs and weight diagnostics
sref  = psw_reference_scale(trial, external)               #    reference scale (pooled outcome SD)
prior = psw_map_prior(s, sigma_ref = sref, tau_c = 16)     # 3. robust PSW-MAP prior
post  = psw_posterior(prior$robust, trial$Y[trial$trt == 0],   # 4. posterior of the treatment effect
                      trial$Y[trial$trt == 1], sigma_ref = sref)
post$effect
post$prob_positive
```

With the same settings and seed, this gives exactly the result of `pswmap(trial, external, tau_c = 16, seed = 1)`.

## Reproducibility

Pass `seed` to `pswmap()` to make an analysis reproducible:

```r
fit = pswmap(trial, external, tau_c = 16, seed = 1)
```

The seed governs every random step: the imputations, the random forest, the XGBoost cross-fitting folds, the MCMC
sampler, the posterior draws and BART. BART is reproducible on one core, which is bartMachine's default. Your
session's random number stream is left unchanged.

The functions of the individual steps (`psw_weights()` and the others) use R's current random number stream, so call
`set.seed()` before them. With BART, also pass `learner_args = list(seed = 1)`, because bartMachine uses its own Java
random number generator.

## Citation

Guo W, Zhu AY. Propensity score integrated MAP prior with historical missingness. Manuscript submitted.

## References

- Neuenschwander B, Capkun-Niggli G, Branson M, Spiegelhalter DJ (2010). Summarizing historical information on
  controls in clinical trials. *Clinical Trials*, 7(1), 5-18.
- Neuenschwander B, Schmidli H (2020). Use of historical data. In: Lesaffre E, Baio G, Boulanger B (eds),
  *Bayesian Methods in Pharmaceutical Research*, 111-137. Chapman and Hall/CRC.
- Neuenschwander B, Weber S, Schmidli H, O'Hagan A (2020). Predictively consistent prior effective sample sizes.
  *Biometrics*, 76(2), 578-587.
- Nguyen TQ, Stuart EA (2024). Multiple imputation for propensity score analysis with covariates missing at random:
  some clarity on "within" and "across" methods. *American Journal of Epidemiology*, 193(10), 1470-1476.
- Rosenbaum PR, Rubin DB (1983). The central role of the propensity score in observational studies for causal
  effects. *Biometrika*, 70(1), 41-55.
- Rubin DB (1976). Inference and missing data. *Biometrika*, 63(3), 581-592.
- Rubin DB (1987). *Multiple Imputation for Nonresponse in Surveys*. Wiley, New York.
- Schmidli H, Gsteiger S, Roychoudhury S, O'Hagan A, Spiegelhalter D, Neuenschwander B (2014). Robust
  meta-analytic-predictive priors in clinical trials with historical control information. *Biometrics*, 70(4),
  1023-1032.
- Weber S, Li Y, Seaman JW, Kakizume T, Schmidli H (2021). Applying meta-analytic-predictive priors with the R
  Bayesian evidence synthesis tools. *Journal of Statistical Software*, 100(19), 1-32.
- White IR, Carlin JB (2010). Bias and efficiency of multiple imputation compared with complete-case analysis for
  missing covariate values. *Statistics in Medicine*, 29(28), 2920-2931.
- Zhu AY, Roy D, Zhu Z, Sailer MO (2024). Propensity score stratified MAP prior and posterior inference for
  incorporating information across multiple potentially heterogeneous data sources. *Journal of Biopharmaceutical
  Statistics*, 34(2), 190-204.

## License

GPL (>= 3). See [LICENSE.md](LICENSE.md).
