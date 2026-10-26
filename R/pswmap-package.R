#' pswmap: propensity-score-weighted MAP priors with missing external covariates
#'
#' @description
#' `pswmap` borrows external control data, such as historical trials or real-world data, for the control arm of a
#' randomized trial. It uses a propensity-score-weighted meta-analytic predictive (PSW-MAP) prior, and it handles
#' baseline covariates that are missing in the external data. The trial itself is assumed to be complete, and the
#' outcome continuous.
#'
#' @section The PSW-MAP analysis:
#' Let \eqn{Z_i = 1} for subjects of the trial and \eqn{Z_i = 0} for subjects of the external sources
#' \eqn{H_1, \dots, H_K}{H_1, ..., H_K}. The analysis has four steps, each available as a function:
#'
#' 1. **Source-specific propensity scores and weights** ([psw_weights()]). For each source \eqn{H_k}, the
#'    propensity score \eqn{e_k(x) = P(Z = 1 \mid X = x)}{e_k(x) = P(Z = 1 | X = x)} is fitted to the trial stacked with \eqn{H_k} only. The
#'    odds weight \eqn{w_i = e_k(x_i) / \{1 - e_k(x_i)\}}{w_i = e_k(x_i) / (1 - e_k(x_i))} reweights \eqn{H_k} to the covariate distribution of the
#'    trial.
#' 2. **Weighted source summaries** ([psw_summaries()]). Each source is summarized by its odds-weighted (Hajek) mean
#'    \eqn{\bar y_k^w}{ybar_k} and design-based variance \eqn{V_k}.
#' 3. **The robust PSW-MAP prior** ([psw_map_prior()]). A normal-normal hierarchical model,
#'    \eqn{\bar y_k^w \sim N(\theta_k, V_k)}{ybar_k ~ N(theta_k, V_k)} with \eqn{\theta_0, \theta_1, \dots, \theta_K \sim N(\mu, \tau^2)}{theta_0, theta_1, ..., theta_K ~ N(mu, tau^2)},
#'    gives the predictive distribution of the trial's control mean \eqn{\theta_0}{theta_0}. A source with poor overlap has
#'    a large \eqn{V_k} and is therefore discounted. The prior is then mixed with a weakly informative component,
#'    which guards against prior-data conflict.
#' 4. **The posterior of the treatment effect** ([psw_posterior()]), \eqn{\beta_1 = \theta_T - \theta_0}{beta_1 = theta_T - theta_0}.
#'
#' [pswmap()] runs all four steps in one call. [psw_reference_scale()] gives the reference scale
#' \eqn{\sigma_{ref}}{sigma_ref} needed when the steps are run one at a time.
#'
#' @section Getting started:
#' Run `pswmap(trial, external, tau_c = 16, seed = 1)`, with the data in the format of [pswmap_example]. The help
#' page of [pswmap()] explains the options, the assumptions and how to read the output. The help pages of the
#' individual steps give the theory behind each step, with formulas and references.
#'
#' @section Missing covariates in the external data:
#' Following Rubin (1976), covariates are missing completely at random (MCAR) when missingness is independent of
#' the data, missing at random (MAR) when it depends only on observed values (other covariates or the outcome), and
#' missing not at random (MNAR) otherwise. Three strategies are available through the argument `missing`:
#' - **Multiple imputation** (`"mi"`, the default): each source is imputed on its own, with the outcome in the
#'   imputation model. Every imputed dataset is analysed up to the weighted source summaries, which are then combined
#'   by Rubin's rules (the within approach). This is valid when the data are missing at random, including when
#'   missingness depends on the outcome.
#' - **Complete-case weighting** (`"cc"`): incomplete external subjects get weight 0. This is consistent when
#'   missingness is independent of the outcome given the covariates, but it discards information.
#' - **Native handling** (`"native"`): the propensity score learner uses the missing values itself, through
#'   missingness indicators, BARTm, or the default split directions of random forests and XGBoost. Because the trial
#'   is complete, a missing value identifies a subject as external. Native methods therefore converge to complete-case
#'   weighting, and they share its bias when missingness depends on the outcome.
#'
#' @section Reproducibility:
#' Pass `seed` to [pswmap()], e.g. `pswmap(trial, external, tau_c = 16, seed = 1)`, to make an analysis
#' reproducible. This covers the imputations, the random forest, the XGBoost cross-fitting folds, the MCMC sampler,
#' the posterior draws and BART (on one core, the bartMachine default). The functions of the individual steps use
#' R's current random number stream, so call [set.seed()] before them. For `learner = "bart"`, also pass
#' `learner_args = list(seed = ...)`, because bartMachine uses its own Java random number generator.
#'
#' @references
#' Guo W, Zhu AY. Propensity score integrated MAP prior with historical missingness. Manuscript submitted.
#'
#' Neuenschwander B, Capkun-Niggli G, Branson M, Spiegelhalter DJ (2010). Summarizing historical information on
#' controls in clinical trials. *Clinical Trials*, 7(1), 5-18.
#'
#' Nguyen TQ, Stuart EA (2024). Multiple imputation for propensity score analysis with covariates missing at random:
#' some clarity on "within" and "across" methods. *American Journal of Epidemiology*, 193(10), 1470-1476.
#'
#' Rubin DB (1976). Inference and missing data. *Biometrika*, 63(3), 581-592.
#'
#' Schmidli H, Gsteiger S, Roychoudhury S, O'Hagan A, Spiegelhalter D, Neuenschwander B (2014). Robust
#' meta-analytic-predictive priors in clinical trials with historical control information. *Biometrics*, 70(4),
#' 1023-1032.
#'
#' Weber S, Li Y, Seaman JW, Kakizume T, Schmidli H (2021). Applying meta-analytic-predictive priors with the R
#' Bayesian evidence synthesis tools. *Journal of Statistical Software*, 100(19), 1-32.
#'
#' @seealso [pswmap()] to run an analysis, and [pswmap_example] for example data.
#' @keywords internal
"_PACKAGE"
