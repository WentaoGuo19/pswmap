#' Source-specific propensity scores and odds weights
#'
#' @description
#' Estimates, for every external source, the propensity score of trial membership and turns it into odds weights
#' that reweight the source to the covariate distribution of the trial. Missing covariates in the external data are
#' handled by multiple imputation, complete-case weighting or the learner's native missing-data handling.
#'
#' This is the first step of [pswmap()]. Most users only need [pswmap()]. Use `psw_weights()` to run and inspect
#' the steps one at a time.
#'
#' @details
#' **Background: weighting external controls to the trial.** External controls are not randomized, so their
#' covariates \eqn{X} follow a different distribution \eqn{f_k(x)} than the trial's \eqn{f_1(x)}. If the outcome
#' depends on \eqn{X} in the same way in both (no unmeasured confounding of source membership), reweighting source
#' \eqn{H_k} by the density ratio \eqn{f_1(x)/f_k(x)} makes its outcome mean an estimate of the control mean of the
#' trial population. The propensity score of trial membership (Rosenbaum and Rubin, 1983) gives this ratio without
#' modelling the densities: in the stacked data of the trial and \eqn{H_k}, Bayes' rule gives
#' \deqn{\frac{e_k(x)}{1 - e_k(x)} = \frac{n_1}{n_k}\, \frac{f_1(x)}{f_k(x)}, \qquad e_k(x) = P(Z = 1 \mid X = x,\ i \in D_1 \cup D_{H_k}),}{e_k(x) / (1 - e_k(x)) = (n_1 / n_k) f_1(x) / f_k(x),   e_k(x) = P(Z = 1 | X = x, subject in the trial or in H_k),}
#' where \eqn{Z = 1} marks trial subjects (both arms) and \eqn{n_1}, \eqn{n_k} are the sample sizes. The constant
#' \eqn{n_1/n_k} cancels in the normalized (Hajek) mean of [psw_summaries()]. This is the "weighting by the odds"
#' used to estimate effects in a target population (Stuart, 2010).
#'
#' **Source-specific scores.** The score is fitted separately for each source, to both arms of the trial stacked
#' with that source only. Each external subject \eqn{i} of \eqn{H_k} gets the odds weight
#' \deqn{w_i = \frac{\min\{\hat e_k(x_i), e_{\max}\}}{1 - \min\{\hat e_k(x_i), e_{\max}\}}.}{w_i = min(e_k(x_i), e_max) / (1 - min(e_k(x_i), e_max)).}
#' The cap \eqn{e_{\max}}{e_max} (`ps_cap`) bounds every weight by \eqn{e_{\max}/(1 - e_{\max})}{e_max/(1 - e_max)}, which is 19 for the
#' default 0.95, so that a few subjects cannot dominate a source. Trial subjects are not weighted. One score fitted
#' to the trial against all sources together (`ps_scope = "pooled"`) has odds proportional to \eqn{f_1/f_H}, where
#' \eqn{f_H} is the density of the pooled external data. It reweights source \eqn{H_k} to \eqn{f_1 f_k / f_H}
#' instead of \eqn{f_1}, so it balances only the mixture of the sources, not each source. The MAP prior would then
#' read the remaining imbalance as between-source heterogeneity.
#'
#' **Background: missing covariates.** Let \eqn{M} be the missingness pattern of a subject (which covariates are
#' missing). Following Rubin (1976), covariates are *missing completely at random* (MCAR) when \eqn{M} is independent
#' of the data, *missing at random* (MAR) when \eqn{M} depends only on observed values (other covariates or the
#' outcome), and *missing not at random* (MNAR) otherwise. Here the trial is complete, and covariates may be missing
#' in the external sources only. The three strategies (`missing`) differ in the assumptions they need:
#' - `"mi"` (default): the missing covariates of each source are multiply imputed with [mice::mice()] (chained
#'   equations, van Buuren and Groothuis-Oudshoorn, 2011), using that source's subjects only and including the
#'   outcome as a predictor. Including the outcome is needed for valid imputations of covariates in an outcome
#'   analysis (Moons et al., 2006). The propensity score is fitted to each of the \eqn{m} completed datasets, giving
#'   one set of weights per imputation, and [psw_summaries()] combines the source summaries by Rubin's rules. This
#'   is the within approach of Nguyen and Stuart (2024). Valid under MAR, including missingness that depends on the
#'   outcome.
#' - `"cc"`: the score is fitted to the trial and the complete external subjects, and incomplete external subjects
#'   get weight 0. Valid when missingness is independent of the outcome given the covariates (White and Carlin,
#'   2010), which covers MCAR and MAR or MNAR that depends on the covariates only. It discards the incomplete
#'   subjects, so it borrows less.
#' - `"native"`: the learner uses the missing values itself, so no imputation model is needed. `"logit"` adds a
#'   missingness indicator for each incomplete covariate, with missing values set to the observed mean. `"bart"`
#'   uses BARTm (Kapelner and Bleich, 2015), and `"rf"` and `"xgb"` send missing values down a learned branch of
#'   each split (missing incorporated in attributes, Twala et al., 2008). Because the trial is complete, a missing
#'   value can only occur in an external subject:
#'   \deqn{P(Z = 1 \mid X_{obs}, M = m) = 0 \quad \text{for every pattern } m \text{ with a missing covariate.}}{P(Z = 1 | X_obs, M = m) = 0 for every pattern m with a missing covariate.}
#'   A learner that estimates this score well gives incomplete subjects weights near 0. Native handling therefore
#'   moves toward complete-case weighting as the missing rate or the sample size grows, and it shares the bias of
#'   complete-case weighting when missingness depends on the outcome.
#'
#' When missingness depends on the outcome, complete-case weighting and native handling shift every source in the
#' same direction. The sources then still agree with one another, so the robust MAP prior cannot detect the bias.
#' Multiple imputation with the outcome in the imputation model avoids it. If the external covariates are complete,
#' `missing` is ignored and the propensity score is fitted to the full data.
#'
#' **Learners** (`learner`). The defaults are those of Guo and Zhu. Override them with `learner_args`.
#' - `"logit"` (default): main-effects logistic regression ([stats::glm()]), fitted probabilities.
#' - `"bart"`: Bayesian additive regression trees (Chipman et al., 2010), a sum of many small regression trees with
#'   a regularizing prior, here a probit model fitted by [bartMachine::bartMachine()] with its defaults (Kapelner
#'   and Bleich, 2016). The score is the in-sample posterior mean probability.
#' - `"rf"`: random forest (Breiman, 2001), an average of trees grown on bootstrap samples, as a probability forest
#'   by [ranger::ranger()] (Wright and Ziegler, 2017) with 500 trees. The score is the out-of-bag probability, the
#'   average over the trees that did not use the subject.
#' - `"xgb"`: gradient-boosted trees (Chen and Guestrin, 2016), trees added one at a time to correct the current
#'   fit, by [xgboost::xgb.train()] with learning rate 0.1, maximum depth 3, minimum child weight 1 and 100 rounds.
#'   The score is obtained by 5-fold cross-fitting stratified by trial membership: each subject's probability comes
#'   from the model fitted to the other folds.
#'
#' In the simulation study of Guo and Zhu, logistic regression gave the lowest mean squared error of the treatment
#' effect. The learners other than `"logit"` need the corresponding package. bartMachine also needs Java, and its
#' memory is set with `options(java.parameters = "-Xmx4g")` before the package is loaded.
#'
#' **Reproducibility.** `psw_weights()` uses R's current random number stream (for the imputations, the random
#' forest and the XGBoost folds), so call [set.seed()] before it for reproducible weights. bartMachine uses its own
#' generator: pass `learner_args = list(seed = ...)`.
#'
#' @param trial Data frame of the randomized trial (both arms), with the outcome, the treatment indicator and the
#'   covariates. The covariates must be numeric (binary covariates coded 0/1, factors expanded into indicator
#'   columns) and complete.
#' @param external The external control data. Either a named list of data frames, one per source, or one data frame
#'   with a column `source` identifying the sources. Each source needs the outcome and the covariates. Covariates may
#'   be missing, but the outcome must be observed.
#' @param outcome Name of the continuous outcome column (default `"Y"`).
#' @param treatment Name of the treatment column of `trial`, coded 1 for treated and 0 for control (default `"trt"`).
#' @param source Name of the column identifying the sources when `external` is one data frame (default `"source"`).
#' @param covariates Character vector of the covariates used in the propensity score. `NULL` (default) uses every
#'   column of `trial` other than the outcome, the treatment and the source.
#' @param missing Handling of missing external covariates: `"mi"` (default), `"cc"` or `"native"`. See Details.
#' @param learner Propensity score learner: `"logit"` (default), `"bart"`, `"rf"` or `"xgb"`. See Details.
#' @param learner_args Named list of arguments that override the learner's settings. These are passed to
#'   [bartMachine::bartMachine()] (e.g. `list(num_trees = 100, seed = 1)`), [ranger::ranger()] (e.g.
#'   `list(num.trees = 1000)`) or [xgboost::xgb.train()] (e.g. `list(eta = 0.05, nrounds = 200)`). For `"xgb"`,
#'   `nfolds` sets the number of cross-fitting folds. Arguments that pswmap sets itself (the data, the outcome, the
#'   missing-data options and the objective) cannot be overridden. Ignored for `"logit"`, with a warning.
#' @param ps_cap Upper cap \eqn{e_{\max}}{e_max} on the propensity score before the weight is formed (default 0.95). `NULL`
#'   applies no cap.
#' @param mi_m Number of imputed datasets \eqn{m} (default 5).
#' @param mi_engine Imputation method of [mice::mice()]: `"default"` (logistic regression for binary covariates and
#'   predictive mean matching for the others), `"cart"` (classification and regression trees) or `"rf"` (random
#'   forests).
#' @param mi_args Named list of further arguments for [mice::mice()], e.g. `list(maxit = 10)`. The arguments `data`,
#'   `m`, `method` and `printFlag` are set by pswmap. The three imputation settings are ignored, with a warning if
#'   they are supplied, when `missing` is `"cc"` or `"native"` or when the external covariates are complete.
#' @param ps_scope `"source"` (default) fits one propensity score per external source. `"pooled"` fits one score for
#'   the trial against all sources together. It does not balance the individual sources and gives a warning. It is
#'   provided for comparison only.
#'
#' @return An object of class `psw_weights`, a list with:
#' \describe{
#'   \item{`sources`}{one element per external source, holding the outcome `y`, the covariates `X` as observed, the
#'     indicator `incomplete` (a missing covariate), and the matrices `ps` (capped propensity scores) and `w` (odds
#'     weights). These matrices have one column per imputed dataset, or a single column without imputation.}
#'   \item{`missing`}{the strategy applied: `"mi"`, `"cc"`, `"native"`, or `"none"` if the external covariates are
#'     complete.}
#'   \item{`m`, `learner`, `ps_scope`, `ps_cap`, `mi_engine`, `covariates`}{the settings used.}
#' }
#' Printing the object shows the weighted source summaries of [psw_summaries()].
#'
#' @references
#' Breiman L (2001). Random forests. *Machine Learning*, 45(1), 5-32.
#'
#' Chen T, Guestrin C (2016). XGBoost: a scalable tree boosting system. In *Proceedings of the 22nd ACM SIGKDD
#' International Conference on Knowledge Discovery and Data Mining*, 785-794.
#'
#' Chipman HA, George EI, McCulloch RE (2010). BART: Bayesian additive regression trees. *The Annals of Applied
#' Statistics*, 4(1), 266-298.
#'
#' Kapelner A, Bleich J (2015). Prediction with missing data via Bayesian additive regression trees. *Canadian
#' Journal of Statistics*, 43(2), 224-239.
#'
#' Kapelner A, Bleich J (2016). bartMachine: machine learning with Bayesian additive regression trees. *Journal of
#' Statistical Software*, 70(4), 1-40.
#'
#' Moons KGM, Donders RART, Stijnen T, Harrell FE (2006). Using the outcome for imputation of missing predictor
#' values was preferred. *Journal of Clinical Epidemiology*, 59(10), 1092-1101.
#'
#' Nguyen TQ, Stuart EA (2024). Multiple imputation for propensity score analysis with covariates missing at random:
#' some clarity on "within" and "across" methods. *American Journal of Epidemiology*, 193(10), 1470-1476.
#'
#' Rosenbaum PR, Rubin DB (1983). The central role of the propensity score in observational studies for causal
#' effects. *Biometrika*, 70(1), 41-55.
#'
#' Rubin DB (1976). Inference and missing data. *Biometrika*, 63(3), 581-592.
#'
#' Stuart EA (2010). Matching methods for causal inference: a review and a look forward. *Statistical Science*,
#' 25(1), 1-21.
#'
#' Twala BETH, Jones MC, Hand DJ (2008). Good methods for coping with missing data in decision trees. *Pattern
#' Recognition Letters*, 29(7), 950-956.
#'
#' van Buuren S, Groothuis-Oudshoorn K (2011). mice: multivariate imputation by chained equations in R. *Journal of
#' Statistical Software*, 45(3), 1-67.
#'
#' White IR, Carlin JB (2010). Bias and efficiency of multiple imputation compared with complete-case analysis for
#' missing covariate values. *Statistics in Medicine*, 29(28), 2920-2931.
#'
#' Wright MN, Ziegler A (2017). ranger: a fast implementation of random forests for high dimensional data in C++
#' and R. *Journal of Statistical Software*, 77(1), 1-17.
#'
#' @seealso [psw_summaries()] for the weighted source summaries, [pswmap()] for the full analysis.
#'
#' @examples
#' data(pswmap_example)
#' w = psw_weights(pswmap_example$trial, pswmap_example$external, missing = "cc")
#' w
#' \donttest{
#' set.seed(1)
#' # multiple imputation with 5 imputed datasets (the default)
#' w_mi = psw_weights(pswmap_example$trial, pswmap_example$external)
#' psw_summaries(w_mi)
#' }
#' @export
psw_weights = function(trial, external, outcome = "Y", treatment = "trt", source = "source", covariates = NULL,
                       missing = c("mi", "cc", "native"), learner = c("logit", "bart", "rf", "xgb"),
                       learner_args = list(), ps_cap = 0.95, mi_m = 5, mi_engine = c("default", "cart", "rf"),
                       mi_args = list(), ps_scope = c("source", "pooled")) {
  if (base::missing(external)) stop("`external` is required: the external control data.", call. = FALSE)
  weights_impl(trial, external, outcome = outcome, treatment = treatment, source = source, covariates = covariates,
               missing = missing, learner = learner, learner_args = learner_args, ps_cap = ps_cap, mi_m = mi_m,
               mi_engine = mi_engine, mi_args = mi_args, ps_scope = ps_scope, given = given_args(match.call()))
}


# Check the settings of the propensity score step (shared by psw_weights() and pswmap(), before any data are used).
# Returns the checked settings.
check_weight_settings = function(missing, learner, learner_args, ps_cap, mi_m, mi_engine, mi_args, ps_scope) {
  missing = check_choice(missing, c("mi", "cc", "native"), "missing")
  learner = check_choice(learner, c("logit", "bart", "rf", "xgb"), "learner")
  mi_engine = check_choice(mi_engine, c("default", "cart", "rf"), "mi_engine")
  ps_scope = check_choice(ps_scope, c("source", "pooled"), "ps_scope")
  ps_cap = check_between(ps_cap, "ps_cap", null_ok = TRUE)
  mi_m = check_whole(mi_m, "mi_m", min = 2)
  learner_args = check_learner_args(check_arg_list(learner_args, "learner_args"), learner)
  mi_args = check_arg_list(mi_args, "mi_args")
  fixed = intersect(names(mi_args), c("data", "m", "method", "printFlag"))
  if (length(fixed))
    stop("`mi_args` must not set ", paste0("`", fixed, "`", collapse = ", "), ": use `mi_m` and `mi_engine`.",
         call. = FALSE)
  list(missing = missing, learner = learner, learner_args = learner_args, ps_cap = ps_cap, mi_m = mi_m,
       mi_engine = mi_engine, mi_args = mi_args, ps_scope = ps_scope)
}


# Check learner_args against the learner: names the learner function accepts, and none of the arguments that the
# package sets itself.
check_learner_args = function(args, learner) {
  if (learner == "logit") return(args)
  pkg = c(bart = "bartMachine", rf = "ranger", xgb = "xgboost")[[learner]]
  need_pkg(pkg, learner)
  if (!length(args)) return(args)
  fixed = switch(learner,
                 bart = c("X", "y", "Xy", "use_missing_data", "use_missing_data_dummies_as_covars", "run_in_sample",
                          "verbose"),
                 rf = c("x", "y", "formula", "data", "probability", "na.action", "dependent.variable.name", "classification"),
                 xgb = c("objective", "data", "label", "params"))
  bad = intersect(names(args), fixed)
  if (length(bad))
    stop("`learner_args` must not set ", paste0("`", bad, "`", collapse = ", "), ": pswmap sets ",
         if (length(bad) == 1) "it" else "them", " for learner = \"", learner, "\".", call. = FALSE)
  known = switch(learner, bart = names(formals(bartMachine::bartMachine)), rf = setdiff(names(formals(ranger::ranger)), "..."),
                 xgb = NULL)
  if (!is.null(known)) {
    unknown = setdiff(names(args), known)
    if (length(unknown))
      stop("unknown `learner_args` for learner = \"", learner, "\": ", paste(unknown, collapse = ", "), " (see ?",
           c(bart = "bartMachine::bartMachine", rf = "ranger::ranger")[[learner]], ").", call. = FALSE)
  }
  if (learner == "xgb") {
    if (!is.null(args$nrounds)) check_whole(args$nrounds, "learner_args$nrounds", min = 1)
    if (!is.null(args$nfolds)) check_whole(args$nfolds, "learner_args$nfolds", min = 2)
  }
  args
}


# The propensity score step. `given`: names of the arguments the caller typed, for warnings about settings that
# have no effect.
weights_impl = function(trial, external, outcome, treatment, source, covariates, missing, learner, learner_args,
                        ps_cap, mi_m, mi_engine, mi_args, ps_scope, given = character(0)) {
  st = check_weight_settings(missing, learner, learner_args, ps_cap, mi_m, mi_engine, mi_args, ps_scope)
  missing = st$missing; learner = st$learner; learner_args = st$learner_args; ps_cap = st$ps_cap; mi_m = st$mi_m
  mi_engine = st$mi_engine; mi_args = st$mi_args; ps_scope = st$ps_scope
  if (learner == "logit" && length(learner_args))
    warning("`learner_args` is ignored for learner = \"logit\" (a main-effects logistic regression).", call. = FALSE)
  if (ps_scope == "pooled")
    warning("ps_scope = \"pooled\" fits one propensity score for the trial against all external sources. It ",
            "balances only their mixture with the trial, not each source, so it is not recommended; ",
            "use the default \"source\".", call. = FALSE)

  dat = prepare_data(trial, external, outcome, treatment, source, covariates, need_external = TRUE)
  cv = dat$covariates
  trial_X = dat$trial[, cv, drop = FALSE]
  ext = dat$external
  ext_X = lapply(ext, function(d) d[, cv, drop = FALSE])
  has_na = any(vapply(ext_X, anyNA, logical(1)))
  strategy = if (!has_na) "none" else missing
  mi_settings = c("mi_m", "mi_engine", "mi_args")
  if (strategy == "none") warn_ignored(c("missing", mi_settings), given, "because the external covariates are complete")
  if (strategy %in% c("cc", "native")) warn_ignored(mi_settings, given, paste0("with missing = \"", strategy, "\""))
  if (strategy == "cc" && any(vapply(ext_X, function(X) sum(stats::complete.cases(X)) < 2, logical(1))))
    stop("complete-case weighting needs at least two complete subjects in every external source.", call. = FALSE)
  if (strategy == "mi" && mi_engine == "rf" && !requireNamespace("ranger", quietly = TRUE))
    stop("mi_engine = \"rf\" needs the package ranger: install.packages(\"ranger\")", call. = FALSE)

  #propensity scores: one column per imputed dataset (MI), one column otherwise
  if (strategy == "mi") {
    imputed = lapply(ext, impute_source, covariates = cv, outcome = dat$outcome, m = mi_m, engine = mi_engine,
                     mi_args = mi_args)
    ps_cols = lapply(seq_len(mi_m), function(i)
      ps_external(trial_X, lapply(imputed, `[[`, i), "full", ps_scope, learner, learner_args))
  } else {
    ps_cols = list(ps_external(trial_X, ext_X, if (strategy == "none") "full" else strategy, ps_scope,
                               learner, learner_args))
  }
  cap = if (is.null(ps_cap)) 1 - 1e-6 else ps_cap   #the PS is always kept below 1
  sources = lapply(names(ext), function(k) {
    ps = do.call(cbind, lapply(ps_cols, `[[`, k))
    ps = pmin(ps, cap)
    list(y = ext[[k]][[dat$outcome]], X = ext_X[[k]], incomplete = !stats::complete.cases(ext_X[[k]]),
         ps = ps, w = ps / (1 - ps))
  })
  names(sources) = names(ext)
  structure(list(sources = sources, missing = strategy, m = length(ps_cols), learner = learner,
                 ps_scope = ps_scope, ps_cap = ps_cap, mi_engine = if (strategy == "mi") mi_engine else NA,
                 covariates = cv, outcome = dat$outcome, treatment = dat$treatment),
            class = "psw_weights")
}


#' @export
print.psw_weights = function(x, ...) {
  cat("Propensity score weights for", length(x$sources), "external source(s)\n")
  cat("  learner:", x$learner, "| scope:", x$ps_scope, "| missing covariates:", x$missing,
      if (x$missing == "mi") paste0("(", x$m, " imputations, engine ", x$mi_engine, ")"),
      "| PS cap:", if (is.null(x$ps_cap)) "none" else x$ps_cap, "\n")
  print(psw_summaries(x), row.names = FALSE, digits = 3)
  invisible(x)
}
