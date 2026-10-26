# Input checking and preparation (internal)

# Split the external data into a named list of sources.
# external: a data frame with a source column, or a list of data frames (one per source).
split_external = function(external, source) {
  if (is.data.frame(external)) {
    if (!source %in% names(external))
      stop("`external` is a data frame, so it needs a column `", source, "` identifying the external sources ",
           "(or pass a list of data frames, one per source).", call. = FALSE)
    ids = external[[source]]
    if (anyNA(ids)) stop("column `", source, "` of `external` has missing values.", call. = FALSE)
    lev = unique(as.character(ids))
    out = lapply(lev, function(s) external[as.character(ids) == s, setdiff(names(external), source), drop = FALSE])
    names(out) = lev
    return(out)
  }
  if (is.list(external) && length(external) > 0 && all(vapply(external, is.data.frame, logical(1)))) {
    if (is.null(names(external)) || any(names(external) == ""))
      names(external) = paste0("H", seq_along(external))
    if (anyDuplicated(names(external))) stop("the external sources must have distinct names.", call. = FALSE)
    return(external)
  }
  stop("`external` must be a data frame with a source column or a list of data frames.", call. = FALSE)
}


# Check the trial and external data and keep only the columns the analysis uses.
# Returns list(trial, external, covariates, outcome, treatment).
prepare_data = function(trial, external, outcome, treatment, source, covariates, need_external) {
  outcome = check_string(outcome, "outcome"); treatment = check_string(treatment, "treatment")
  source = check_string(source, "source")
  if (outcome == treatment) stop("`outcome` and `treatment` must be different columns.", call. = FALSE)
  if (!is.data.frame(trial)) stop("`trial` must be a data frame.", call. = FALSE)
  for (v in c(outcome, treatment))
    if (!v %in% names(trial)) stop("column `", v, "` is not in `trial`.", call. = FALSE)

  #outcome: continuous only
  y = trial[[outcome]]
  if (!is.numeric(y)) stop("the outcome `", outcome, "` must be numeric (continuous).", call. = FALSE)
  if (anyNA(y)) stop("the outcome has missing values in `trial`.", call. = FALSE)
  if (any(!is.finite(y))) stop("the outcome has infinite values in `trial`.", call. = FALSE)
  if (length(unique(y)) <= 2)
    stop("pswmap supports continuous outcomes only; `", outcome, "` has at most two distinct values.", call. = FALSE)

  #treatment: 0/1 with both arms present
  tr = trial[[treatment]]
  if (is.logical(tr)) tr = as.numeric(tr)
  if (!is.numeric(tr) || anyNA(tr) || !all(tr %in% c(0, 1)) || !all(c(0, 1) %in% tr))
    stop("the treatment `", treatment, "` must be coded 0 (control) / 1 (treated), with both arms present.", call. = FALSE)
  for (a in 0:1) {
    arm = if (a == 0) "control" else "treated"
    if (sum(tr == a) < 2) stop("the ", arm, " arm of `trial` needs at least two subjects.", call. = FALSE)
    if (stats::var(y[tr == a]) == 0) stop("the outcomes of the ", arm, " arm of `trial` are all equal.", call. = FALSE)
  }
  trial[[treatment]] = tr

  #external sources
  ext = NULL
  if (!is.null(external)) ext = split_external(external, source)
  if (!is.null(ext) && length(ext) == 0) stop("`external` contains no subjects.", call. = FALSE)
  if (need_external && is.null(ext)) stop("borrowing needs `external` data.", call. = FALSE)

  #covariates: by default all trial columns except outcome, treatment and source
  if (is.null(covariates)) covariates = setdiff(names(trial), c(outcome, treatment, source))
  if (!is.character(covariates) || anyNA(covariates))
    stop("`covariates` must be NULL or a character vector of column names.", call. = FALSE)
  if (anyDuplicated(covariates))
    stop("`covariates` contains duplicates: ", paste(unique(covariates[duplicated(covariates)]), collapse = ", "),
         call. = FALSE)
  bad = intersect(covariates, c(outcome, treatment))
  if (length(bad)) stop("`covariates` must not contain the outcome or the treatment: ", paste(bad, collapse = ", "),
                        call. = FALSE)
  if (need_external && length(covariates) == 0) stop("no covariates for the propensity score.", call. = FALSE)
  miss_trial = setdiff(covariates, names(trial))
  if (length(miss_trial)) stop("covariates not in `trial`: ", paste(miss_trial, collapse = ", "), call. = FALSE)
  for (v in covariates) {
    if (is.logical(trial[[v]])) trial[[v]] = as.numeric(trial[[v]])
    if (!is.numeric(trial[[v]]))
      stop("covariate `", v, "` must be numeric (code binary covariates as 0/1 and expand factors into ",
           "indicator columns).", call. = FALSE)
    if (any(is.infinite(trial[[v]]))) stop("covariate `", v, "` has infinite values in `trial`.", call. = FALSE)
  }
  if (need_external && anyNA(trial[, covariates]))
    stop("the trial has missing covariates; pswmap assumes a complete trial (missing values are allowed in the ",
         "external data only).", call. = FALSE)
  trial = trial[, c(outcome, treatment, covariates), drop = FALSE]

  if (!is.null(ext)) {
    ext = lapply(names(ext), function(k) {
      d = ext[[k]]
      need = c(outcome, covariates)
      miss = setdiff(need, names(d))
      if (length(miss)) stop("external source `", k, "` lacks columns: ", paste(miss, collapse = ", "), call. = FALSE)
      if (!is.numeric(d[[outcome]]) || anyNA(d[[outcome]]) || any(!is.finite(d[[outcome]])))
        stop("the outcome must be numeric and observed (finite) for every subject of external source `", k, "`.",
             call. = FALSE)
      if (nrow(d) < 2) stop("external source `", k, "` has fewer than two subjects.", call. = FALSE)
      if (stats::var(d[[outcome]]) == 0)
        stop("the outcomes of external source `", k, "` are all equal.", call. = FALSE)
      for (v in covariates) {
        if (is.logical(d[[v]])) d[[v]] = as.numeric(d[[v]])
        if (!is.numeric(d[[v]])) stop("covariate `", v, "` of external source `", k, "` must be numeric.", call. = FALSE)
        if (any(is.infinite(d[[v]])))
          stop("covariate `", v, "` has infinite values in external source `", k, "`.", call. = FALSE)
        if (need_external && all(is.na(d[[v]])))
          stop("covariate `", v, "` is missing for every subject of external source `", k, "`; drop it from ",
               "`covariates` or from the analysis.", call. = FALSE)
      }
      d[, need, drop = FALSE]
    })
    names(ext) = names(split_external(external, source))
  }
  list(trial = trial, external = ext, covariates = covariates, outcome = outcome, treatment = treatment)
}


# Reference scale: pooled (degrees-of-freedom-weighted) outcome SD of the trial's control arm and every external
# source, computed from the unweighted outcomes so that it does not depend on the weights.
pooled_sd = function(trial_control_y, external_y) {
  groups = c(list(trial_control_y), external_y)
  df = vapply(groups, function(y) length(y) - 1, numeric(1))
  v = vapply(groups, stats::var, numeric(1))
  sqrt(sum(df * v) / sum(df))
}


#' Reference scale of the PSW-MAP analysis
#'
#' @description
#' Computes the reference scale \eqn{\sigma_{ref}}{sigma_ref}: the pooled outcome standard deviation of the trial's
#' control arm and every external source. It sets the scale of the hyperpriors, of the robust component and of the
#' treated-arm prior. [pswmap()] computes it automatically. Use this function when running the steps one at a time
#' with [psw_map_prior()] and [psw_posterior()].
#'
#' @details
#' With groups \eqn{g = 0, 1, \dots, K}{g = 0, 1, ..., K} (the trial's control arm and the \eqn{K} sources), group
#' sizes \eqn{n_g} and unweighted outcome variances \eqn{s_g^2},
#' \deqn{\sigma_{ref}^2 = \frac{\sum_{g=0}^K (n_g - 1)\, s_g^2}{\sum_{g=0}^K (n_g - 1)}.}{sigma_ref^2 = sum_g (n_g - 1) s_g^2 / sum_g (n_g - 1).}
#' Because it uses the unweighted outcomes, it does not depend on the propensity score weights.
#'
#' The reference scale is the standard deviation of a single outcome, and it makes the settings of the analysis
#' free of the outcome's units (Neuenschwander and Schmidli, 2020). The prior on the between-source standard
#' deviation \eqn{\tau}{tau} is set relative to it (scale \eqn{\sigma_{ref}/c}{sigma_ref/c}). A normal component with standard deviation
#' \eqn{\sigma_{ref}}{sigma_ref} carries the information of one subject, which is used for the robust component and the treated-arm
#' prior. The prior effective sample size counts subjects with this standard deviation.
#'
#' @references
#' Neuenschwander B, Schmidli H (2020). Use of historical data. In: Lesaffre E, Baio G, Boulanger B (eds),
#' *Bayesian Methods in Pharmaceutical Research*, 111-137. Chapman and Hall/CRC.
#'
#' @inheritParams psw_weights
#' @param external The external control data, as in [psw_weights()], or `NULL` to use the trial's control arm only.
#'
#' @return The reference scale, a positive number.
#' @seealso [psw_map_prior()], [psw_posterior()]
#' @examples
#' data(pswmap_example)
#' psw_reference_scale(pswmap_example$trial, pswmap_example$external)
#' @export
psw_reference_scale = function(trial, external = NULL, outcome = "Y", treatment = "trt", source = "source") {
  dat = prepare_data(trial, external, outcome, treatment, source, covariates = character(0), need_external = FALSE)
  y = dat$trial[[outcome]]
  pooled_sd(y[dat$trial[[treatment]] == 0],
            if (is.null(dat$external)) list() else lapply(dat$external, `[[`, outcome))
}
