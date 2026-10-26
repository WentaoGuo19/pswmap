#' Borrow external controls through a propensity-score-weighted MAP prior
#'
#' @description
#' `pswmap()` analyses a randomized trial whose control arm is augmented with one or more external control sources.
#' It fits source-specific propensity scores, summarizes each source by its odds-weighted mean, builds a robust
#' meta-analytic predictive (MAP) prior for the trial's control mean, and returns the posterior of the treatment
#' effect. Covariates may be missing in the external data, but not in the trial. With `borrow = FALSE` the trial
#' is analysed on its own, which gives the reference analysis without borrowing.
#'
#' This is the function most users need. Advanced users who want to compute and inspect each step themselves can
#' run the steps one at a time with [psw_weights()], [psw_summaries()], [psw_reference_scale()], [psw_map_prior()]
#' and [psw_posterior()]. With the same settings and seed, this gives exactly the same result.
#'
#' @details
#' **The analysis.** For each external source \eqn{H_k}, the propensity score \eqn{e_k(x) = P(Z = 1 \mid X = x)}{e_k(x) = P(Z = 1 | X = x)} is
#' fitted to the trial (\eqn{Z = 1}) stacked with \eqn{H_k} (\eqn{Z = 0}). Each external subject gets the odds
#' weight \eqn{w_i = e_k(x_i)/\{1 - e_k(x_i)\}}{w_i = e_k(x_i)/(1 - e_k(x_i))}, with the score capped at `ps_cap`. Source \eqn{H_k} is summarized
#' by its weighted mean \eqn{\bar y_k^w}{ybar_k} and design-based variance \eqn{V_k} (see [psw_summaries()]). The PSW-MAP
#' prior is the predictive distribution of the control mean \eqn{\theta_0}{theta_0} in the hierarchical model
#' \deqn{\bar y_k^w \mid \theta_k \sim N(\theta_k, V_k), \qquad \theta_0, \theta_1, \dots, \theta_K \mid \mu, \tau
#'   \sim N(\mu, \tau^2),}{ybar_k | theta_k ~ N(theta_k, V_k),   theta_0, theta_1, ..., theta_K | mu, tau ~ N(mu, tau^2),}
#' with \eqn{\tau \sim}{tau ~} half-N\eqn{(0, (\sigma_{ref}/c)^2)}{(0, (sigma_ref/c)^2)}. It is mixed with a weakly informative component of
#' weight `robust_weight` (see [psw_map_prior()]), updated with the trial's control arm, and combined with the
#' treated arm's posterior to give the posterior of \eqn{\beta_1 = \theta_T - \theta_0}{beta_1 = theta_T - theta_0} (see [psw_posterior()]).
#'
#' **Assumptions.** The analysis relies on three assumptions, each addressed by one part of the method:
#' - *Measured covariates explain the differences.* Given the covariates, the control outcome of an external
#'   subject has the same distribution as in the trial. The propensity score weights adjust for the measured
#'   differences. Remaining differences between sources (for example in standard of care) are absorbed by the
#'   between-source standard deviation \eqn{\tau}{tau} and, under prior-data conflict, by the robust component.
#' - *Overlap.* Every source contains subjects resembling the trial subjects. With poor overlap the weights are
#'   very unequal, the source's variance \eqn{V_k} is large, and the source is discounted. The Kish effective
#'   sample size in `sources` shows the overlap of each source.
#' - *A missing-data assumption that matches the chosen strategy* (see below).
#'
#' **Missing external covariates** (`missing`). Covariates may be missing completely at random (MCAR), at random
#' given the observed data (MAR), or not at random (MNAR). [psw_weights()] gives the definitions and details.
#' - `"mi"` (recommended) imputes each source on its own, with the outcome in the imputation model. It analyses
#'   every imputed dataset up to the source summaries and combines these by Rubin's rules. Valid under MAR,
#'   including when missingness depends on the outcome.
#' - `"cc"` gives incomplete external subjects weight 0. Valid when missingness is independent of the outcome given
#'   the covariates, at the cost of less borrowing.
#' - `"native"` lets the learner handle the missing values itself. Because the trial is complete, a missing value
#'   identifies a subject as external, so these methods move toward complete-case weighting as the missing rate or
#'   the sample size grows, and share its assumptions.
#'
#' Use `"mi"` unless you have reason to believe that missingness does not depend on the outcome. Then `"cc"` is a
#' simple alternative that needs no imputation model. The relative weight of incomplete subjects in `sources`
#' shows how much weight incomplete subjects keep under each strategy.
#'
#' **Choosing `tau_c`.** The heterogeneity constant \eqn{c} sets the scale \eqn{\sigma_{ref}/c}{sigma_ref/c} of the half-normal
#' prior on the between-source standard deviation \eqn{\tau}{tau}. Values 32, 16, 8 and 4 correspond to small,
#' moderate, substantial and large anticipated heterogeneity (Zhu et al., 2024, following Neuenschwander and
#' Schmidli, 2020). A smaller \eqn{c} borrows less and protects better against prior-data conflict. There is no
#' default, because the choice should reflect what is known about the sources.
#'
#' **Reading the output.** The prior effective sample size (ESS) states how many trial control subjects the external
#' data are worth. Compare the posterior with that of `borrow = FALSE`, the trial alone: the gain from borrowing
#' shows as a narrower credible interval. A much smaller ESS than the external sample sizes signals poor overlap or heterogeneity
#' between the sources.
#'
#' **Defaults.** All other settings default to the values recommended by Guo and Zhu:
#' - a main-effects logistic propensity score, capped at 0.95;
#' - multiple imputation with 5 imputed datasets;
#' - a robust component of weight 0.2 with the information of one subject;
#' - the MCMC and mixture settings of [map_control()].
#'
#' **Reproducibility.** Set `seed` to make an analysis reproducible. The seed then governs every random step:
#' - the imputations;
#' - the random forest and the XGBoost cross-fitting folds;
#' - the MCMC sampler and the posterior draws;
#' - BART, through bartMachine's own random number generator.
#'
#' Your session's random number stream is restored afterwards. With the default `seed = NULL` the analysis continues
#' the current stream, so calling [set.seed()] beforehand works as well, except for BART. bartMachine is
#' reproducible only on one core, its default.
#'
#' @inheritParams psw_weights
#' @param borrow `TRUE` (default) borrows from the external data. `FALSE` analyses the trial alone: both arms get
#'   the same vague prior \eqn{N(\bar y, (100\,\sigma_{ref})^2)}{N(ybar, (100 sigma_ref)^2)}, centred at the trial mean \eqn{\bar y}{ybar}, and
#'   `external` may be omitted. If given, `external` is used only for the default `sigma_ref`. The settings of the
#'   borrowing steps (`covariates`, `tau_c`, the propensity score and imputation settings, `map`, the robust
#'   component and `treated_prior`) are then ignored, with a warning if they are supplied.
#' @param tau_c Heterogeneity constant \eqn{c} of the half-normal prior on \eqn{\tau}{tau}, whose scale is
#'   \eqn{\sigma_{ref}/c}{sigma_ref/c}. Required when `borrow = TRUE`. See Details.
#' @param sigma_ref Reference scale \eqn{\sigma_{ref}}{sigma_ref}, which sets the scales of the hyperpriors, the robust
#'   component and the treated-arm prior. `NULL` (default) uses the pooled outcome standard deviation of the trial's
#'   control arm and every external source,
#'   \eqn{\sigma_{ref}^2 = \sum_g (n_g - 1) s_g^2 / \sum_g (n_g - 1)}{sigma_ref^2 = sum_g (n_g - 1) s_g^2 / sum_g (n_g - 1)}. This is computed from the unweighted outcomes,
#'   so it does not depend on the weights.
#' @param map MCMC and mixture settings of the MAP prior, from [map_control()].
#' @param robust_weight,robust_mean,robust_sd Weight, mean and standard deviation of the robust component. See
#'   [psw_map_prior()].
#' @param treated_prior Prior of the treated-arm mean, an RBesT normal mixture (see [RBesT::mixnorm()]). `NULL`
#'   (default) uses the weakly informative \eqn{N(m, \sigma_{ref}^2)}{N(m, sigma_ref^2)}, where \eqn{m} is the mean of the robust
#'   prior.
#' @param n_draws,level,threshold Number of posterior draws per arm, credible level, and decision threshold for
#'   \eqn{P(\beta_1 > 0 \mid \text{data})}{P(beta_1 > 0 | data)}. See [psw_posterior()].
#' @param seed Optional single number. If given, the analysis is run with this random seed and is reproducible,
#'   including BART (unless `learner_args` sets bartMachine's own `seed`). The caller's random number stream is
#'   left unchanged. `NULL` (default) uses the current stream.
#'
#' @return An object of class `pswmap`, a list with:
#' \describe{
#'   \item{`effect`}{posterior mean, SD, median and credible limits of the treatment effect \eqn{\beta_1}{beta_1}.}
#'   \item{`prob_positive`, `effective`}{\eqn{P(\beta_1 > 0 \mid \text{data})}{P(beta_1 > 0 | data)} and whether it exceeds `threshold`.}
#'   \item{`sources`}{the weighted source summaries from [psw_summaries()].}
#'   \item{`prior`}{the MAP prior, the robust prior, the prior effective sample size and the MCMC fit, from
#'     [psw_map_prior()].}
#'   \item{`weights`}{the `psw_weights` object.}
#'   \item{`post_control`, `post_treated`}{the posterior mixtures of the two arm means.}
#'   \item{`sigma_ref`, `borrow`, `call`}{the reference scale, whether external data were borrowed, and the call.}
#' }
#' With `borrow = FALSE`, `sources`, `prior` and `weights` are `NULL`. Use `print()` and `summary()` to display the
#' results.
#'
#' @references
#' Guo W, Zhu AY. Propensity score integrated MAP prior with historical missingness. Manuscript submitted.
#'
#' Neuenschwander B, Schmidli H (2020). Use of historical data. In: Lesaffre E, Baio G, Boulanger B (eds),
#' *Bayesian Methods in Pharmaceutical Research*, 111-137. Chapman and Hall/CRC.
#'
#' Schmidli H, Gsteiger S, Roychoudhury S, O'Hagan A, Spiegelhalter D, Neuenschwander B (2014). Robust
#' meta-analytic-predictive priors in clinical trials with historical control information. *Biometrics*, 70(4),
#' 1023-1032.
#'
#' Zhu AY, Roy D, Zhu Z, Sailer MO (2024). Propensity score stratified MAP prior and posterior inference for
#' incorporating information across multiple potentially heterogeneous data sources. *Journal of Biopharmaceutical
#' Statistics*, 34(2), 190-204.
#'
#' @seealso [psw_weights()], [psw_summaries()], [psw_map_prior()], [psw_posterior()], [map_control()],
#'   [pswmap_example]
#'
#' @examples
#' \donttest{
#' data(pswmap_example)
#' trial = pswmap_example$trial          # 100 treated, 50 control
#' external = pswmap_example$external    # three external sources with missing covariates
#'
#' # PSW-MAP with multiple imputation (default), anticipating moderate heterogeneity
#' fit = pswmap(trial, external, tau_c = 16, seed = 1)
#' fit
#' summary(fit)
#'
#' # complete-case weighting, and the trial alone for reference
#' fit_cc = pswmap(trial, external, tau_c = 16, missing = "cc", seed = 1)
#' fit_trial = pswmap(trial, borrow = FALSE, seed = 1)
#'
#' # sensitivity to the anticipated heterogeneity: with the same seed, the weights are the same
#' fit_large = pswmap(trial, external, tau_c = 4, seed = 1)
#' }
#' @export
pswmap = function(trial, external = NULL, outcome = "Y", treatment = "trt", source = "source", covariates = NULL,
                  borrow = TRUE, tau_c, missing = c("mi", "cc", "native"),
                  learner = c("logit", "bart", "rf", "xgb"), learner_args = list(), ps_cap = 0.95,
                  mi_m = 5, mi_engine = c("default", "cart", "rf"), mi_args = list(), ps_scope = c("source", "pooled"),
                  sigma_ref = NULL, map = map_control(), robust_weight = 0.2, robust_mean = NULL, robust_sd = NULL,
                  treated_prior = NULL, n_draws = 10000, level = 0.95, threshold = 0.975,
                  seed = NULL) {
  cl = match.call()
  given = given_args(cl)
  #check every setting before any computation
  borrow = check_flag(borrow, "borrow")
  if (!is.null(seed)) check_number(seed, "seed")
  check_positive(sigma_ref, "sigma_ref", null_ok = TRUE)
  n_draws = check_whole(n_draws, "n_draws", min = 2)
  check_between(level, "level"); check_between(threshold, "threshold")
  borrow_only = c("covariates", "tau_c", "missing", "learner", "learner_args", "ps_cap", "mi_m", "mi_engine",
                  "mi_args", "ps_scope", "map", "robust_weight", "robust_mean", "robust_sd", "treated_prior")
  if (borrow) {
    if (base::missing(tau_c))
      stop("`tau_c` is required when borrow = TRUE: the heterogeneity constant c of the prior on tau ",
           "(scale sigma_ref / c), e.g. 32, 16, 8 or 4 for small, moderate, substantial or large heterogeneity.",
           call. = FALSE)
    check_positive(tau_c, "tau_c")
    st = check_weight_settings(missing, learner, learner_args, ps_cap, mi_m, mi_engine, mi_args, ps_scope)
    learner = st$learner; learner_args = st$learner_args
    if (!inherits(map, "psw_map_control")) stop("`map` must come from map_control().", call. = FALSE)
    check_robust(robust_weight, robust_mean, robust_sd, warn = FALSE)   #psw_map_prior() warns
    check_normmix(treated_prior, "treated_prior", null_ok = TRUE)
  } else {
    warn_ignored(borrow_only, given, "when borrow = FALSE (the trial is analysed alone)")
  }

  if (!is.null(seed)) {
    #set the seed for this analysis, and restore the caller's random number stream on exit
    old_seed = get0(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(restore_seed(old_seed), add = TRUE)
    set.seed(seed)
    #bartMachine draws its random numbers in Java: pass the seed to it as well, unless the user gave one
    if (borrow && learner == "bart" && is.null(learner_args$seed)) learner_args$seed = seed
  }
  dat = prepare_data(trial, external, outcome, treatment, source, covariates, need_external = borrow)
  y = dat$trial[[outcome]]; tr = dat$trial[[treatment]]
  if (is.null(sigma_ref))
    sigma_ref = psw_reference_scale(trial, external, outcome = outcome, treatment = treatment, source = source)

  if (!borrow) {
    vague = RBesT::mixnorm(c(1, mean(y), 100 * sigma_ref), sigma = sigma_ref)
    post = psw_posterior(vague, y[tr == 0], y[tr == 1], sigma_ref, treated_prior = vague, n_draws = n_draws,
                         level = level, threshold = threshold)
    return(structure(c(post, list(borrow = FALSE, sources = NULL, prior = NULL, weights = NULL,
                                  sigma_ref = sigma_ref, call = cl)), class = "pswmap"))
  }

  weights = weights_impl(trial, external, outcome = outcome, treatment = treatment, source = source,
                         covariates = dat$covariates, missing = missing, learner = learner,
                         learner_args = learner_args, ps_cap = ps_cap, mi_m = mi_m, mi_engine = mi_engine,
                         mi_args = mi_args, ps_scope = ps_scope, given = given)
  sources = psw_summaries(weights)
  prior = psw_map_prior(sources, sigma_ref = sigma_ref, tau_c = tau_c, control = map,
                        robust_weight = robust_weight, robust_mean = robust_mean, robust_sd = robust_sd)
  post = psw_posterior(prior$robust, y[tr == 0], y[tr == 1], sigma_ref, treated_prior = treated_prior,
                       n_draws = n_draws, level = level, threshold = threshold)
  structure(c(post, list(borrow = TRUE, sources = sources, prior = prior, weights = weights, sigma_ref = sigma_ref,
                         call = cl)), class = "pswmap")
}


#' @export
print.pswmap = function(x, digits = 3, ...) {
  if (x$borrow) {
    w = x$weights
    cat("PSW-MAP analysis:", length(w$sources), "external source(s), learner", w$learner,
        "| missing covariates:", w$missing,
        if (w$missing == "mi") paste0("(", w$m, " imputations)"), "\n")
    cat("Heterogeneity constant c =", x$prior$tau_c, "| reference scale =", format(x$sigma_ref, digits = digits),
        "| prior ESS =", format(x$prior$ess, digits = digits), "\n\nExternal sources:\n")
    print(x$sources, row.names = FALSE, digits = digits)
  } else {
    cat("Analysis of the trial alone (no borrowing)\n")
  }
  e = x$effect
  cat("\nTreatment effect (treated - control):\n")
  cat("  posterior mean ", format(e["mean"], digits = digits), " (SD ", format(e["sd"], digits = digits), ") ",
      100 * x$level, "% CrI [", format(e["lower"], digits = digits), ", ", format(e["upper"], digits = digits), "]",
      sep = "")
  cat("\n  P(effect > 0) =", format(x$prob_positive, digits = digits),
      "->", if (x$effective) "effective" else "not effective", "at threshold", x$threshold, "\n")
  invisible(x)
}


#' @export
summary.pswmap = function(object, ...) {
  out = list(effect = object$effect, prob_positive = object$prob_positive, effective = object$effective,
             sources = object$sources, prior_ess = if (object$borrow) object$prior$ess else 0,
             map_prior = if (object$borrow) object$prior$map else NULL,
             robust_prior = if (object$borrow) object$prior$robust else NULL)
  class(out) = "summary.pswmap"
  out
}


#' @export
print.summary.pswmap = function(x, digits = 3, ...) {
  print(round(x$effect, digits))
  cat("P(effect > 0) =", format(x$prob_positive, digits = digits), "| prior ESS =", format(x$prior_ess, digits = digits), "\n")
  if (!is.null(x$sources)) { cat("\n"); print(x$sources, row.names = FALSE, digits = digits) }
  if (!is.null(x$robust_prior)) { cat("\nRobust PSW-MAP prior:\n"); print(x$robust_prior) }
  invisible(x)
}

# Restore the caller's random number stream (or its absence) after an analysis run with `seed`.
restore_seed = function(old) {
  if (is.null(old)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  } else assign(".Random.seed", old, envir = globalenv())
}

utils::globalVariables(c("Y_bar_w", "Y_bar_w.se", "Study", "n"))
