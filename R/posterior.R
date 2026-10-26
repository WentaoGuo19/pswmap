#' Posterior of the treatment effect
#'
#' @description
#' Updates the priors of the two arm means with the trial data and summarizes the treatment effect
#' \eqn{\beta_1 = \theta_T - \theta_0}{beta_1 = theta_T - theta_0}, the difference between the treated and control means.
#'
#' This is the last step of [pswmap()]. Most users only need [pswmap()].
#'
#' @details
#' **Conjugate updating of a mixture prior.** Each arm's data enter through the sample mean \eqn{\bar y}{ybar} and its
#' standard error \eqn{s_e = s/\sqrt{n}}{s_e = s / sqrt(n)}, with \eqn{\bar y \mid \theta \sim N(\theta, s_e^2)}{ybar | theta ~ N(theta, s_e^2)}. A normal-mixture prior
#' \eqn{\sum_j \pi_j N(m_j, s_j^2)}{sum_j pi_j N(m_j, s_j^2)} then gives a normal-mixture posterior with components
#' \deqn{s_j'^2 = \left(\frac{1}{s_j^2} + \frac{1}{s_e^2}\right)^{-1}, \qquad m_j' = s_j'^2 \left(\frac{m_j}{s_j^2} +
#'   \frac{\bar y}{s_e^2}\right), \qquad \pi_j' \propto \pi_j\, N(\bar y;\, m_j,\, s_j^2 + s_e^2).}{s_j'^2 = 1 / (1/s_j^2 + 1/s_e^2),   m_j' = s_j'^2 (m_j/s_j^2 + ybar/s_e^2),   pi_j' proportional to pi_j N(ybar; m_j, s_j^2 + s_e^2).}
#' The new weights \eqn{\pi_j'}{pi_j'} favour the components that predict the observed mean well. For the robust prior of
#' [psw_map_prior()], this is how a conflict between the trial and the external sources moves weight to the
#' weakly informative component (dynamic borrowing).
#'
#' **Treatment effect.** The two arms are updated independently ([RBesT::postmix()]), and the posterior of
#' \eqn{\beta_1}{beta_1} is summarized by `n_draws` independent draws from each posterior. No external data inform the treated
#' arm. Its default prior is the weakly informative \eqn{N(m, \sigma_{ref}^2)}{N(m, sigma_ref^2)}, where \eqn{m} is the mean of the
#' control-arm prior, which carries the information of one subject.
#'
#' **Decision rule.** The treatment is declared effective when \eqn{P(\beta_1 > 0 \mid \text{data}) > }{P(beta_1 > 0 | data) >} `threshold`.
#' The default threshold 0.975 corresponds to the lower limit of the two-sided 95% credible interval lying above 0,
#' the Bayesian analogue of a one-sided test at level 0.025. It assumes that larger outcomes are better. For
#' smaller-is-better outcomes, multiply the outcome by -1 before the analysis.
#'
#' @param control_prior Prior of the control mean \eqn{\theta_0}{theta_0}, an RBesT normal mixture, for example the `robust`
#'   element of [psw_map_prior()].
#' @param y_control,y_treated Outcomes of the trial's control and treated arms.
#' @param sigma_ref Reference scale \eqn{\sigma_{ref}}{sigma_ref}, used by the default treated-arm prior.
#' @param treated_prior Prior of the treated mean \eqn{\theta_T}{theta_T}, an RBesT normal mixture (see
#'   [RBesT::mixnorm()]). `NULL` (default) uses \eqn{N(m, \sigma_{ref}^2)}{N(m, sigma_ref^2)}. See Details.
#' @param n_draws Number of posterior draws per arm (default 10,000).
#' @param level Credible level of the interval (default 0.95).
#' @param threshold Decision threshold for \eqn{P(\beta_1 > 0 \mid \text{data})}{P(beta_1 > 0 | data)} (default 0.975).
#'
#' @return A list with:
#' \describe{
#'   \item{`effect`}{posterior mean, SD, median, and lower and upper credible limits of \eqn{\beta_1}{beta_1}.}
#'   \item{`prob_positive`}{\eqn{P(\beta_1 > 0 \mid \text{data})}{P(beta_1 > 0 | data)}.}
#'   \item{`effective`}{`TRUE` if `prob_positive` exceeds `threshold`.}
#'   \item{`post_control`, `post_treated`}{the posterior mixtures of \eqn{\theta_0}{theta_0} and \eqn{\theta_T}{theta_T}.}
#'   \item{`level`, `threshold`}{the settings used.}
#' }
#'
#' @references
#' Schmidli H, Gsteiger S, Roychoudhury S, O'Hagan A, Spiegelhalter D, Neuenschwander B (2014). Robust
#' meta-analytic-predictive priors in clinical trials with historical control information. *Biometrics*, 70(4),
#' 1023-1032.
#'
#' Weber S, Li Y, Seaman JW, Kakizume T, Schmidli H (2021). Applying meta-analytic-predictive priors with the R
#' Bayesian evidence synthesis tools. *Journal of Statistical Software*, 100(19), 1-32.
#'
#' @seealso [psw_map_prior()], [pswmap()]
#'
#' @examples
#' \donttest{
#' data(pswmap_example)
#' tr = pswmap_example$trial
#' s = psw_summaries(psw_weights(tr, pswmap_example$external, missing = "cc"))
#' set.seed(1)
#' sref = psw_reference_scale(pswmap_example$trial, pswmap_example$external)
#' prior = psw_map_prior(s, sigma_ref = sref, tau_c = 16)
#' post = psw_posterior(prior$robust, tr$Y[tr$trt == 0], tr$Y[tr$trt == 1], sigma_ref = sref)
#' post$effect
#' post$prob_positive
#' }
#' @export
psw_posterior = function(control_prior, y_control, y_treated, sigma_ref, treated_prior = NULL, n_draws = 10000,
                         level = 0.95, threshold = 0.975) {
  check_normmix(control_prior, "control_prior")
  check_arm(y_control, "y_control"); check_arm(y_treated, "y_treated")
  check_positive(sigma_ref, "sigma_ref")
  check_normmix(treated_prior, "treated_prior", null_ok = TRUE)
  n_draws = check_whole(n_draws, "n_draws", min = 2)
  check_between(level, "level"); check_between(threshold, "threshold")
  if (is.null(treated_prior))
    treated_prior = RBesT::mixnorm(c(1, sum(control_prior[1, ] * control_prior[2, ]), sigma_ref), sigma = sigma_ref)
  post_c = suppressMessages(RBesT::postmix(control_prior, data = y_control))
  post_t = suppressMessages(RBesT::postmix(treated_prior, data = y_treated))
  draws = RBesT::rmix(post_t, n = n_draws) - RBesT::rmix(post_c, n = n_draws)
  a = (1 - level) / 2
  q = stats::quantile(draws, c(a, 0.5, 1 - a), names = FALSE)
  effect = c(mean = mean(draws), sd = stats::sd(draws), median = q[2], lower = q[1], upper = q[3])
  p = mean(draws > 0)
  list(effect = effect, prob_positive = p, effective = p > threshold, post_control = post_c, post_treated = post_t,
       level = level, threshold = threshold)
}
