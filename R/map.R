#' Settings of the MAP prior
#'
#' @description
#' Collects the settings of the meta-analytic predictive (MAP) prior: the prior on the overall mean \eqn{\mu}{mu}, the
#' MCMC settings of [RBesT::gMAP()], and the normal-mixture approximation of the prior. The heterogeneity constant
#' `tau_c` is not set here. It is an argument of [pswmap()] and [psw_map_prior()]. The defaults are those used by
#' Guo and Zhu.
#'
#' @details
#' **MCMC.** The hierarchical model of [psw_map_prior()] has no closed-form posterior, so [RBesT::gMAP()] samples it
#' with the Hamiltonian Monte Carlo sampler of Stan. Each chain runs `iter` iterations, of which the first `warmup`
#' tune the sampler and are discarded, and every `thin`-th draw of the rest is kept. With few sources, the posterior
#' of \eqn{\tau}{tau} has much mass near 0, where the sampler can produce divergent transitions. A target acceptance
#' rate `adapt_delta` close to 1 makes the sampler take smaller steps and avoids them.
#'
#' **Mixture approximation.** The MCMC draws of \eqn{\theta_0}{theta_0} are approximated by a mixture of normal
#' densities, \eqn{\sum_j \pi_j N(m_j, s_j^2)}{sum_j pi_j N(m_j, s_j^2)}, fitted by the EM algorithm for each number of components in
#' `mix_components` ([RBesT::automixfit()]). The fit with the lowest penalized Akaike information criterion is kept.
#' A mixture is a conjugate prior for a normal mean, so the posterior in [psw_posterior()] is again a mixture with
#' closed-form weights, means and variances.
#'
#' @param mu_mean,mu_sd Mean and standard deviation of the normal prior on the overall mean,
#'   \eqn{\mu \sim N(\mu_0, \sigma_\mu^2)}{mu ~ N(mu_0, sigma_mu^2)}. The default mean is 0. `mu_sd = NULL` (default) uses
#'   \eqn{\sigma_\mu = 100\,\sigma_{ref}}{sigma_mu = 100 sigma_ref}, which is effectively flat over any plausible range of \eqn{\theta_0}{theta_0}.
#' @param iter,warmup,thin,chains MCMC settings of [RBesT::gMAP()]. The defaults (4 chains of 6,000 iterations, of
#'   which 2,000 are warm-up, thinned by 4) give 4,000 posterior draws. These are the RBesT defaults.
#' @param adapt_delta Target acceptance rate of the Stan sampler (default 0.99). A high value reduces divergent
#'   transitions, which are common in hierarchical models with few sources.
#' @param cores Number of cores used to run the chains in parallel (default `getOption("mc.cores", 1)`).
#' @param mix_components Numbers of components tried when the MCMC sample of \eqn{\theta_0}{theta_0} is approximated by a
#'   normal mixture with [RBesT::automixfit()] (default 1 to 4). The fit with the best penalized AIC is kept.
#' @param min_sd_ratio Candidate mixtures containing a component whose standard deviation is below `min_sd_ratio`
#'   times the standard deviation of the whole mixture are excluded (default 0.05). The EM algorithm occasionally
#'   collapses a component onto a single extreme MCMC draw. Such a component leaves the prior essentially unchanged
#'   but inflates its effective sample size.
#'
#' @return A list of class `psw_map_control`, to be passed as `map` to [pswmap()] or as `control` to
#'   [psw_map_prior()].
#'
#' @seealso [psw_map_prior()]
#'
#' @examples
#' map_control()                                  # the defaults
#' map_control(iter = 10000, warmup = 5000, cores = 4)
#' @export
map_control = function(mu_mean = 0, mu_sd = NULL, iter = 6000, warmup = 2000, thin = 4, chains = 4,
                       adapt_delta = 0.99, cores = getOption("mc.cores", 1L), mix_components = 1:4,
                       min_sd_ratio = 0.05) {
  mu_mean = check_number(mu_mean, "mu_mean")
  mu_sd = check_positive(mu_sd, "mu_sd", null_ok = TRUE)
  iter = check_whole(iter, "iter", min = 2)
  warmup = check_whole(warmup, "warmup", min = 1)
  if (warmup >= iter) stop("`warmup` must be smaller than `iter`.", call. = FALSE)
  thin = check_whole(thin, "thin", min = 1)
  chains = check_whole(chains, "chains", min = 1)
  if ((iter - warmup) * chains / thin < 100)
    stop("too few MCMC draws: (iter - warmup) * chains / thin must be at least 100.", call. = FALSE)
  adapt_delta = check_between(adapt_delta, "adapt_delta")
  cores = check_whole(cores, "cores", min = 1)
  if (!(is.numeric(mix_components) && length(mix_components) >= 1 && !anyNA(mix_components) &&
        all(mix_components == round(mix_components)) && all(mix_components >= 1) && !anyDuplicated(mix_components)))
    stop("`mix_components` must be distinct whole numbers of at least 1, e.g. 1:4.", call. = FALSE)
  min_sd_ratio = check_between(min_sd_ratio, "min_sd_ratio", closed_lower = TRUE)
  structure(list(mu_mean = mu_mean, mu_sd = mu_sd, iter = iter, warmup = warmup, thin = thin, chains = chains,
                 adapt_delta = adapt_delta, cores = cores, mix_components = mix_components,
                 min_sd_ratio = min_sd_ratio), class = "psw_map_control")
}


#' The robust propensity-score-weighted MAP prior
#'
#' @description
#' Builds the prior for the trial's control mean from the weighted source summaries: a meta-analytic predictive
#' (MAP) prior, made robust by a weakly informative component.
#'
#' This is the third step of [pswmap()]. Most users only need [pswmap()].
#'
#' @details
#' **Background: the MAP approach.** The MAP approach (Neuenschwander et al., 2010) treats the control means of the
#' external sources and of the trial as exchangeable: similar, but not identical. A hierarchical model estimates how
#' much they vary, through the between-source standard deviation \eqn{\tau}{tau}, and predicts the control mean of
#' the trial from the sources. The amount of borrowing follows from the data: sources that agree give a small
#' \eqn{\tau}{tau} and an informative prior, and sources that disagree give a large \eqn{\tau}{tau} and a vague one.
#' Unlike a power prior, which discounts pooled external data by one fixed parameter, the MAP prior accounts for
#' several sources and their heterogeneity directly.
#'
#' **The PSW-MAP model.** The weighted means \eqn{\bar y_k^w}{ybar_k} and variances \eqn{V_k} of [psw_summaries()] enter the
#' normal-normal hierarchical model
#' \deqn{\bar y_k^w \mid \theta_k \sim N(\theta_k, V_k), \qquad \theta_0, \theta_1, \dots, \theta_K \mid \mu, \tau
#'   \overset{iid}{\sim} N(\mu, \tau^2),}{ybar_k | theta_k ~ N(theta_k, V_k),   theta_0, theta_1, ..., theta_K | mu, tau ~ iid N(mu, tau^2),}
#' with \eqn{\mu \sim N(\mu_0, \sigma_\mu^2)}{mu ~ N(mu_0, sigma_mu^2)} and \eqn{\tau \sim}{tau ~} half-N\eqn{(0, (\sigma_{ref}/c)^2)}{(0, (sigma_ref/c)^2)}. Here
#' \eqn{\theta_0}{theta_0} is the control mean of the trial and \eqn{\theta_k}{theta_k} the covariate-adjusted control mean implied
#' by source \eqn{H_k}. Marginally, \eqn{\bar y_k^w \sim N(\mu, \tau^2 + V_k)}{ybar_k ~ N(mu, tau^2 + V_k)}, so a source with poor overlap (large
#' \eqn{V_k}) has less influence. The PSW-MAP prior is the posterior predictive distribution of \eqn{\theta_0}{theta_0},
#' \deqn{p(\theta_0) = \int\!\!\int p(\theta_0 \mid \mu, \tau)\, p(\mu, \tau \mid \bar y_1^w, \dots, \bar y_K^w)\,
#'   d\mu\, d\tau.}{p(theta_0) = integral of p(theta_0 | mu, tau) p(mu, tau | ybar_1, ..., ybar_K) d(mu) d(tau).}
#' It is obtained by MCMC with [RBesT::gMAP()] (Weber et al., 2021) and approximated by a normal mixture with
#' [RBesT::automixfit()] (see [map_control()]).
#'
#' **The heterogeneity constant `tau_c`.** With few sources, the data say little about \eqn{\tau}{tau}, so its prior
#' matters. The half-normal prior has scale \eqn{\sigma_{ref}/c}{sigma_ref/c}, which expresses the anticipated
#' between-source standard deviation relative to the outcome standard deviation. Values \eqn{c} = 32, 16, 8 and 4
#' correspond to small, moderate, substantial and large anticipated heterogeneity (Zhu et al., 2024, adapting
#' Neuenschwander and Schmidli, 2020). A smaller \eqn{c} allows more heterogeneity, borrows less, and protects
#' better against prior-data conflict. Choose it before seeing the trial data, and report the sensitivity of the
#' conclusion to it.
#'
#' **The robust prior.** If the trial's control mean differs from what the sources predict (prior-data conflict),
#' an informative prior biases the estimate. The robust prior of Schmidli et al. (2014) adds a weakly informative
#' component:
#' \deqn{p_{robust}(\theta_0) = (1 - w_R)\, p(\theta_0) + w_R\, N(\theta_0; m_R, s_R^2).}{p_robust(theta_0) = (1 - w_R) p(theta_0) + w_R N(theta_0; m_R, s_R^2).}
#' Under conflict, the trial data shift the posterior weight toward this component, so the external information is
#' discounted automatically. By default, \eqn{w_R = 0.2}, \eqn{m_R} is the mean of the MAP prior, and
#' \eqn{s_R = \sigma_{ref}}{s_R = sigma_ref}, which carries the information of a single subject.
#'
#' **Prior effective sample size.** The ESS states how many trial control subjects the prior is worth. It is
#' computed by the expected local-information-ratio method (Neuenschwander et al., 2020),
#' \deqn{\text{ESS} = E_{p(\theta_0)}\left\{-\frac{\partial^2 \log p(\theta_0)}{\partial \theta_0^2}\right\}\, \sigma_{ref}^2,}{ESS = E over p(theta_0) of { -d^2 log p(theta_0) / d theta_0^2 } * sigma_ref^2,}
#' the expected curvature of the log prior relative to the information \eqn{1/\sigma_{ref}^2}{1/sigma_ref^2} of one subject. For a
#' normal prior \eqn{N(m, s^2)} this is \eqn{\sigma_{ref}^2/s^2}{sigma_ref^2/s^2}.
#'
#' @param summaries Data frame from [psw_summaries()].
#' @param sigma_ref Reference scale \eqn{\sigma_{ref}}{sigma_ref}, e.g. the pooled outcome standard deviation of the trial's
#'   control arm and the external sources (the default of [pswmap()]).
#' @param tau_c Heterogeneity constant \eqn{c}: the half-normal prior on \eqn{\tau}{tau} has scale \eqn{\sigma_{ref}/c}{sigma_ref/c}.
#'   No default. Values 32, 16, 8 and 4 correspond to small, moderate, substantial and large anticipated
#'   heterogeneity (Zhu et al., 2024). See Details.
#' @param control MCMC and mixture settings from [map_control()].
#' @param robust_weight Weight \eqn{w_R} of the robust component (default 0.2). With 0 the MAP prior is returned
#'   without a robust component.
#' @param robust_mean Mean \eqn{m_R} of the robust component. `NULL` (default) uses the mean of the MAP prior.
#' @param robust_sd Standard deviation \eqn{s_R} of the robust component. `NULL` (default) uses `sigma_ref`.
#'
#' @return A list with:
#' \describe{
#'   \item{`map`}{the PSW-MAP prior, an RBesT normal mixture with reference scale `sigma_ref`.}
#'   \item{`robust`}{the robust PSW-MAP prior.}
#'   \item{`ess`}{the prior effective sample size of the robust prior (expected local-information-ratio method).}
#'   \item{`gmap`}{the [RBesT::gMAP()] fit, for MCMC diagnostics.}
#'   \item{`tau_c`, `sigma_ref`}{the settings used.}
#' }
#'
#' @references
#' Neuenschwander B, Capkun-Niggli G, Branson M, Spiegelhalter DJ (2010). Summarizing historical information on
#' controls in clinical trials. *Clinical Trials*, 7(1), 5-18.
#'
#' Neuenschwander B, Schmidli H (2020). Use of historical data. In: Lesaffre E, Baio G, Boulanger B (eds),
#' *Bayesian Methods in Pharmaceutical Research*, 111-137. Chapman and Hall/CRC.
#'
#' Neuenschwander B, Weber S, Schmidli H, O'Hagan A (2020). Predictively consistent prior effective sample sizes.
#' *Biometrics*, 76(2), 578-587.
#'
#' Schmidli H, Gsteiger S, Roychoudhury S, O'Hagan A, Spiegelhalter D, Neuenschwander B (2014). Robust
#' meta-analytic-predictive priors in clinical trials with historical control information. *Biometrics*, 70(4),
#' 1023-1032.
#'
#' Weber S, Li Y, Seaman JW, Kakizume T, Schmidli H (2021). Applying meta-analytic-predictive priors with the R
#' Bayesian evidence synthesis tools. *Journal of Statistical Software*, 100(19), 1-32.
#'
#' Zhu AY, Roy D, Zhu Z, Sailer MO (2024). Propensity score stratified MAP prior and posterior inference for
#' incorporating information across multiple potentially heterogeneous data sources. *Journal of Biopharmaceutical
#' Statistics*, 34(2), 190-204.
#'
#' @seealso [psw_summaries()], [map_control()], [psw_posterior()]
#'
#' @examples
#' \donttest{
#' data(pswmap_example)
#' s = psw_summaries(psw_weights(pswmap_example$trial, pswmap_example$external, missing = "cc"))
#' set.seed(1)
#' sref = psw_reference_scale(pswmap_example$trial, pswmap_example$external)
#' prior = psw_map_prior(s, sigma_ref = sref, tau_c = 16)
#' prior$robust
#' prior$ess
#' }
#' @export
psw_map_prior = function(summaries, sigma_ref, tau_c, control = map_control(), robust_weight = 0.2,
                         robust_mean = NULL, robust_sd = NULL) {
  if (missing(tau_c)) stop("`tau_c` is required: the heterogeneity constant c of the prior on tau (scale sigma_ref / c), ",
                           "e.g. 32, 16, 8 or 4 for small, moderate, substantial or large heterogeneity.", call. = FALSE)
  check_positive(tau_c, "tau_c")
  check_positive(sigma_ref, "sigma_ref")
  if (!inherits(control, "psw_map_control")) stop("`control` must come from map_control().", call. = FALSE)
  check_robust(robust_weight, robust_mean, robust_sd)
  summaries = check_summaries(summaries)
  mu_sd = if (is.null(control$mu_sd)) 100 * sigma_ref else control$mu_sd

  gdat = data.frame(Study = seq_len(nrow(summaries)), n = summaries$kish_ess, Y_bar_w = summaries$mean,
                    Y_bar_w.se = summaries$se)
  old = options(RBesT.MC.control = list(adapt_delta = control$adapt_delta))
  on.exit(options(old), add = TRUE)
  fit = suppressWarnings(RBesT::gMAP(cbind(Y_bar_w, Y_bar_w.se) ~ 1 | Study, data = gdat, weights = n,
                                     family = stats::gaussian, beta.prior = cbind(control$mu_mean, mu_sd),
                                     tau.dist = "HalfNormal", tau.prior = cbind(0, sigma_ref / tau_c),
                                     iter = control$iter, warmup = control$warmup, chains = control$chains,
                                     thin = control$thin, cores = control$cores))
  map = nondegenerate_mix(RBesT::automixfit(fit, Nc = control$mix_components), control$min_sd_ratio)
  #reference scale of the mixture: sigma_ref, which does not depend on the weights (gMAP estimates its own scale
  #from the weighted summaries). It sets the robust component, the effective sample size and the treated-arm prior
  RBesT::sigma(map) = sigma_ref

  robust = map
  if (robust_weight > 0) {
    rs = if (is.null(robust_sd)) sigma_ref else robust_sd
    robust = suppressMessages(if (is.null(robust_mean)) RBesT::robustify(map, weight = robust_weight, n = (sigma_ref / rs)^2) else
      RBesT::robustify(map, weight = robust_weight, mean = robust_mean, n = (sigma_ref / rs)^2))
  }
  list(map = map, robust = robust, ess = suppressMessages(RBesT::ess(robust)), gmap = fit, tau_c = tau_c,
       sigma_ref = sigma_ref)
}


# Check the robust component; with warn = TRUE, warn when its mean or SD is set without a robust component.
check_robust = function(robust_weight, robust_mean, robust_sd, warn = TRUE) {
  check_between(robust_weight, "robust_weight", closed_lower = TRUE)
  if (!is.null(robust_mean)) check_number(robust_mean, "robust_mean")
  check_positive(robust_sd, "robust_sd", null_ok = TRUE)
  set = c("robust_mean", "robust_sd")[c(!is.null(robust_mean), !is.null(robust_sd))]
  if (warn && robust_weight == 0) warn_ignored(set, set, "because robust_weight = 0")
}


# Check the source summaries: a data frame with the columns mean, se and kish_ess (as from psw_summaries()).
check_summaries = function(summaries) {
  if (!is.data.frame(summaries) || nrow(summaries) == 0)
    stop("`summaries` must be a data frame with one row per external source, e.g. from psw_summaries().", call. = FALSE)
  miss = setdiff(c("mean", "se", "kish_ess"), names(summaries))
  if (length(miss)) stop("`summaries` lacks columns: ", paste(miss, collapse = ", "), " (see psw_summaries()).",
                         call. = FALSE)
  for (v in c("mean", "se", "kish_ess"))
    if (!is.numeric(summaries[[v]]) || anyNA(summaries[[v]]) || any(!is.finite(summaries[[v]])))
      stop("column `", v, "` of `summaries` must be numeric and finite.", call. = FALSE)
  if (any(summaries$se <= 0)) stop("the standard errors `se` in `summaries` must be positive.", call. = FALSE)
  if (any(summaries$kish_ess <= 0)) stop("`kish_ess` in `summaries` must be positive.", call. = FALSE)
  if (nrow(summaries) == 1)
    warning("only one external source: the between-source standard deviation tau is informed by its prior alone, ",
            "so the amount of borrowing rests on `tau_c`.", call. = FALSE)
  summaries
}


# Exclude degenerate normal-mixture fits: automixfit() keeps the best fit by penalized AIC, but EM occasionally
# collapses a component onto one extreme MCMC draw. Return the best candidate whose components all have a standard
# deviation of at least min_sd_ratio times that of the whole mixture (the one-component fit always qualifies).
nondegenerate_mix = function(fit, min_sd_ratio = 0.05) {
  ok = function(m) {
    w = m[1, ]; mu = m[2, ]; s = m[3, ]
    sd_mix = sqrt(sum(w * (s^2 + mu^2)) - sum(w * mu)^2)
    all(s >= min_sd_ratio * sd_mix)
  }
  if (ok(fit)) return(fit)
  for (m in attr(fit, "models")) if (ok(m)) return(m)
  stop("no candidate mixture without a degenerate component.", call. = FALSE)
}
