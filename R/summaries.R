#' Weighted summaries of the external sources
#'
#' @description
#' Summarizes every external source by its odds-weighted mean and design-based variance. These are the data of the
#' meta-analytic predictive model in [psw_map_prior()]. Also reports two diagnostics of the weights.
#'
#' This is the second step of [pswmap()]. Most users only need [pswmap()].
#'
#' @details
#' **Weighted mean and its variance.** For source \eqn{H_k} with outcomes \eqn{y_i} and odds weights \eqn{w_i}
#' from [psw_weights()], the weighted mean is the Hajek (normalized) estimator of the control mean of the trial
#' population, and its variance is the linearization (sandwich) variance of this ratio estimator with the weights
#' treated as fixed:
#' \deqn{\bar y_k^w = \frac{\sum_{i \in H_k} w_i y_i}{\sum_{i \in H_k} w_i}, \qquad
#'   V_k = \frac{\sum_{i \in H_k} w_i^2 (y_i - \bar y_k^w)^2}{(\sum_{i \in H_k} w_i)^2}.}{ybar_k = sum(w_i y_i) / sum(w_i),   V_k = sum(w_i^2 (y_i - ybar_k)^2) / (sum(w_i))^2.}
#' Because the weights are normalized, their overall scale cancels. A source whose subjects need very unequal
#' weights to resemble the trial (poor overlap) has a large \eqn{V_k}, so the MAP prior automatically gives it less
#' influence.
#'
#' **Multiple imputation.** With \eqn{m} imputed datasets, the summaries \eqn{(\bar y_k^{w(j)}, V_k^{(j)})}{(ybar_k(j), V_k(j))},
#' \eqn{j = 1, \dots, m}{j = 1, ..., m}, are combined by Rubin's rules (Rubin, 1987):
#' \deqn{\bar y_k^w = \frac{1}{m}\sum_{j=1}^m \bar y_k^{w(j)}, \qquad V_k = \bar W_k + \left(1 + \frac{1}{m}\right)
#'   B_k,}{ybar_k = mean over j of ybar_k(j),   V_k = W_k + (1 + 1/m) B_k,}
#' where \eqn{\bar W_k}{W_k} is the mean of the within-imputation variances \eqn{V_k^{(j)}} and \eqn{B_k} the variance of the
#' \eqn{\bar y_k^{w(j)}}{ybar_k(j)} between imputations. The term \eqn{B_k} carries the uncertainty about the missing covariates
#' into the amount of borrowing, and the factor \eqn{1 + 1/m} accounts for the finite number of imputations.
#'
#' **Diagnostics.** The Kish effective sample size
#' \deqn{\text{ESS}_k = \frac{(\sum_i w_i)^2}{\sum_i w_i^2}}{ESS_k = (sum w_i)^2 / sum w_i^2}
#' is the number of equally weighted subjects that would give the same precision (Kish, 1965). A value far below
#' the source size signals poor overlap with the trial. The relative weight of incomplete subjects is their share
#' of the total weight divided by their share of the subjects,
#' \deqn{\text{RW}_k = \frac{\sum_{i \in I_k} w_i / \sum_{i \in H_k} w_i}{|I_k| / n_k},}{RW_k = (sum of w_i over I_k / sum of w_i over H_k) / (|I_k| / n_k),}
#' where \eqn{I_k} is the set of subjects with at least one missing covariate. It is 0 under complete-case
#' weighting and shows how far native handling has moved toward it. When missingness is completely at random, the
#' correct value is about 1. When missingness depends on the covariates, the correct value differs from 1, because
#' the incomplete subjects then differ in their covariates. Both diagnostics are averaged over the imputations.
#'
#' @param weights An object of class `psw_weights` from [psw_weights()].
#'
#' @return A data frame with one row per external source and the columns:
#' \describe{
#'   \item{`source`}{the source name.}
#'   \item{`n`, `n_incomplete`}{the number of subjects, and of subjects with at least one missing covariate.}
#'   \item{`mean`, `se`}{the weighted mean \eqn{\bar y_k^w}{ybar_k} and its standard error \eqn{\sqrt{V_k}}{sqrt(V_k)}.}
#'   \item{`kish_ess`}{the Kish effective sample size of the weights (averaged over imputations).}
#'   \item{`relweight_incomplete`}{the relative weight of incomplete subjects (averaged over imputations), `NA` if
#'     the source is complete.}
#' }
#'
#' @references
#' Kish L (1965). *Survey Sampling*. Wiley, New York.
#'
#' Rubin DB (1987). *Multiple Imputation for Nonresponse in Surveys*. Wiley, New York.
#'
#' @seealso [psw_weights()], [psw_map_prior()]
#'
#' @examples
#' data(pswmap_example)
#' w = psw_weights(pswmap_example$trial, pswmap_example$external, missing = "cc")
#' psw_summaries(w)
#' @export
psw_summaries = function(weights) {
  if (!inherits(weights, "psw_weights")) stop("`weights` must come from psw_weights().", call. = FALSE)
  rows = lapply(names(weights$sources), function(k) {
    s = weights$sources[[k]]
    y = s$y
    per = apply(s$w, 2, function(w) {
      ybar = sum(w * y) / sum(w)
      c(mean = ybar, var = sum(w^2 * (y - ybar)^2) / sum(w)^2, kish = sum(w)^2 / sum(w^2),
        relw = if (any(s$incomplete)) (sum(w[s$incomplete]) / sum(w)) / mean(s$incomplete) else NA_real_)
    })
    m = ncol(per)
    v = if (m == 1) per["var", 1] else mean(per["var", ]) + (1 + 1 / m) * stats::var(per["mean", ])   #Rubin
    data.frame(source = k, n = length(y), n_incomplete = sum(s$incomplete), mean = mean(per["mean", ]),
               se = sqrt(v), kish_ess = mean(per["kish", ]), relweight_incomplete = mean(per["relw", ]))
  })
  do.call(rbind, rows)
}
