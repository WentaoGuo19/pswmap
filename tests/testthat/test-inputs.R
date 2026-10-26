test_that("tau_c is required when borrowing", {
  expect_error(pswmap(ex_trial(), ex_ext()), "tau_c")
})

test_that("only continuous outcomes are accepted", {
  tr = ex_trial(); tr$Y = as.numeric(tr$Y > 1)
  expect_error(pswmap(tr, ex_ext(), tau_c = 16), "continuous")
})

test_that("the trial must be complete and the treatment coded 0/1", {
  tr = ex_trial(); tr$X1[1] = NA
  expect_error(psw_weights(tr, ex_ext()), "complete trial")
  tr = ex_trial(); tr$trt = tr$trt + 1
  expect_error(psw_weights(tr, ex_ext()), "0 \\(control\\)")
})

test_that("covariates must be numeric and present in every source", {
  tr = ex_trial(); tr$X1 = factor(tr$X1)
  expect_error(psw_weights(tr, ex_ext()), "numeric")
  ext = ex_ext(); ext$H2$X3 = NULL
  expect_error(psw_weights(ex_trial(), ext), "lacks columns")
})

test_that("external data can be a list or one data frame with a source column", {
  ext = ex_ext()
  df = do.call(rbind, lapply(names(ext), function(k) cbind(ext[[k]], source = k)))
  a = psw_summaries(psw_weights(ex_trial(), ext, missing = "cc"))
  b = psw_summaries(psw_weights(ex_trial(), df, missing = "cc"))
  expect_equal(a, b)
  expect_error(psw_weights(ex_trial(), df[, names(df) != "source"]), "source")
  expect_equal(psw_summaries(psw_weights(ex_trial(), unname(ext), missing = "cc"))$source, c("H1", "H2", "H3"))
})

test_that("the reference scale is the pooled SD of the trial controls and the sources", {
  tr = ex_trial(); ext = ex_ext()
  groups = c(list(tr$Y[tr$trt == 0]), lapply(ext, `[[`, "Y"))
  df = sapply(groups, length) - 1
  expect_equal(psw_reference_scale(tr, ext), sqrt(sum(df * sapply(groups, var)) / sum(df)))
  expect_equal(psw_reference_scale(tr), sd(tr$Y[tr$trt == 0]))
})

test_that("column names and choices are checked", {
  expect_error(psw_weights(ex_trial(), ex_ext(), outcome = c("Y", "X1")), "single column name")
  expect_error(psw_weights(ex_trial(), ex_ext(), outcome = "trt"), "different columns")
  expect_error(psw_weights(ex_trial(), ex_ext(), missing = "impute"), "`missing` must be one of")
  expect_error(psw_weights(ex_trial(), ex_ext(), learner = "glm"), "`learner` must be one of")
  expect_error(psw_weights(ex_trial(), ex_ext(), mi_engine = "pmm"), "`mi_engine` must be one of")
  expect_error(psw_weights(ex_trial(), ex_ext(), ps_scope = 2), "`ps_scope` must be one of")
  expect_error(psw_weights(ex_trial(), ex_ext(), covariates = c("X1", "X1")), "duplicates")
  expect_error(psw_weights(ex_trial(), ex_ext(), covariates = c("X1", "Y")), "must not contain")
  expect_error(psw_weights(ex_trial(), ex_ext(), covariates = 1:3), "character vector")
  expect_error(psw_weights(ex_trial()), "`external` is required")
})

test_that("the data are checked", {
  tr = ex_trial(); tr$Y[1] = Inf
  expect_error(psw_weights(tr, ex_ext()), "infinite")
  tr = ex_trial(); tr$trt = c(1, rep(0, nrow(tr) - 1))
  expect_error(psw_weights(tr, ex_ext()), "treated arm of `trial` needs at least two")
  ext = ex_ext(); ext$H1$Y = 1
  expect_error(psw_weights(ex_trial(), ext), "all equal")
  ext = ex_ext(); ext$H2$X2 = NA_real_
  expect_error(psw_weights(ex_trial(), ext), "missing for every subject")
  ext = ex_ext(); ext$H3$X5[1] = -Inf
  expect_error(psw_weights(ex_trial(), ext), "infinite")
  expect_error(psw_weights(ex_trial(), list()), "external")
})

test_that("propensity score settings are checked", {
  expect_error(psw_weights(ex_trial(), ex_ext(), ps_cap = 1), "ps_cap")
  expect_error(psw_weights(ex_trial(), ex_ext(), ps_cap = NA), "ps_cap")
  expect_error(psw_weights(ex_trial(), ex_ext(), mi_m = 1), "mi_m")
  expect_error(psw_weights(ex_trial(), ex_ext(), mi_m = 2.5), "mi_m")
  expect_error(psw_weights(ex_trial(), ex_ext(), learner_args = 3), "must be a list")
  expect_error(psw_weights(ex_trial(), ex_ext(), mi_args = list(10)), "distinct name")
  expect_error(psw_weights(ex_trial(), ex_ext(), mi_args = list(m = 10)), "mi_m")
  skip_if_not_installed("ranger")
  expect_error(psw_weights(ex_trial(), ex_ext(), learner = "rf", learner_args = list(num.tree = 10)), "unknown")
  expect_error(psw_weights(ex_trial(), ex_ext(), learner = "rf", learner_args = list(probability = FALSE)),
               "must not set")
  skip_if_not_installed("xgboost")
  expect_error(psw_weights(ex_trial(), ex_ext(), learner = "xgb", learner_args = list(nfolds = 1)), "nfolds")
})

test_that("settings without effect give a warning", {
  expect_warning(psw_weights(ex_trial(), ex_ext(), learner_args = list(maxit = 5), missing = "cc"),
                 "ignored for learner = \"logit\"")
  expect_warning(psw_weights(ex_trial(), ex_ext(), missing = "cc", mi_m = 10), "`mi_m` is ignored")
  ext = lapply(ex_ext(), function(d) d[stats::complete.cases(d), ])
  expect_warning(psw_weights(ex_trial(), ext, missing = "cc"), "covariates are complete")
  expect_no_warning(psw_weights(ex_trial(), ext))
  expect_warning(pswmap(ex_trial(), borrow = FALSE, tau_c = 16, learner = "rf"),
                 "`tau_c`, `learner` are ignored when borrow = FALSE")
  expect_no_warning(pswmap(ex_trial(), ex_ext(), borrow = FALSE))
})

test_that("pswmap() checks its settings before any computation", {
  expect_error(pswmap(ex_trial(), ex_ext(), tau_c = 16, borrow = "yes"), "TRUE or FALSE")
  expect_error(pswmap(ex_trial(), ex_ext(), tau_c = -1), "tau_c")
  expect_error(pswmap(ex_trial(), ex_ext(), tau_c = 16, sigma_ref = 0), "sigma_ref")
  expect_error(pswmap(ex_trial(), ex_ext(), tau_c = 16, map = list(iter = 100)), "map_control")
  expect_error(pswmap(ex_trial(), ex_ext(), tau_c = 16, robust_weight = 1), "robust_weight")
  expect_error(pswmap(ex_trial(), ex_ext(), tau_c = 16, robust_sd = -2), "robust_sd")
  expect_error(pswmap(ex_trial(), ex_ext(), tau_c = 16, treated_prior = c(0, 1)), "normal mixture")
  expect_error(pswmap(ex_trial(), ex_ext(), tau_c = 16, level = 95), "level")
  expect_error(pswmap(ex_trial(), ex_ext(), tau_c = 16, n_draws = 0), "n_draws")
  expect_error(pswmap(ex_trial(), ex_ext(), tau_c = 16, seed = "a"), "seed")
  expect_error(pswmap(ex_trial(), ex_ext(), tau_c = 16, missing = "x"), "`missing` must be one of")
})

test_that("map_control() checks its settings", {
  expect_error(map_control(iter = 1000, warmup = 1000), "smaller than `iter`")
  expect_error(map_control(chains = 0), "chains")
  expect_error(map_control(thin = 1.5), "thin")
  expect_error(map_control(adapt_delta = 1), "adapt_delta")
  expect_error(map_control(mu_sd = 0), "mu_sd")
  expect_error(map_control(mu_mean = NA), "mu_mean")
  expect_error(map_control(mix_components = c(0, 1)), "mix_components")
  expect_error(map_control(min_sd_ratio = 1), "min_sd_ratio")
  expect_error(map_control(iter = 60, warmup = 10, chains = 1, thin = 1), "too few MCMC draws")
})

test_that("psw_map_prior() and psw_posterior() check their inputs", {
  s = data.frame(source = "H1", mean = 2, se = 0.1, kish_ess = 50)
  expect_error(psw_map_prior(s[, -2], sigma_ref = 1, tau_c = 16), "lacks columns")
  expect_error(psw_map_prior(transform(s, se = 0), sigma_ref = 1, tau_c = 16), "positive")
  expect_error(psw_map_prior(list(), sigma_ref = 1, tau_c = 16), "data frame")
  expect_error(psw_map_prior(s, sigma_ref = 1, tau_c = 16, robust_mean = "a"), "robust_mean")
  expect_error(psw_map_prior(s, sigma_ref = 1, tau_c = 16, control = list()), "map_control")
  p = RBesT::mixnorm(c(1, 0, 1), sigma = 1)
  y = c(0.1, 0.5, 0.9)
  expect_error(psw_posterior(c(0, 1), y, y, sigma_ref = 1), "control_prior")
  expect_error(psw_posterior(p, 1, y, sigma_ref = 1), "y_control")
  expect_error(psw_posterior(p, y, c(y, NA), sigma_ref = 1), "y_treated")
  expect_error(psw_posterior(p, y, rep(1, 3), sigma_ref = 1), "all equal")
  expect_error(psw_posterior(p, y, y, sigma_ref = -1), "sigma_ref")
  expect_error(psw_posterior(p, y, y, sigma_ref = 1, threshold = 1), "threshold")
  expect_error(psw_posterior(p, y, y, sigma_ref = 1, treated_prior = RBesT::mixbeta(c(1, 1, 1))), "treated_prior")
  expect_s3_class(psw_posterior(p, y, y, sigma_ref = 1)$post_control, "normMix")
})

test_that("a single source and robust settings without a robust component give warnings", {
  s = data.frame(source = "H1", mean = 2, se = 0.1, kish_ess = 50)
  expect_warning(expect_warning(
    psw_map_prior(s, sigma_ref = 1, tau_c = 16, control = fast_map(), robust_weight = 0, robust_sd = 2),
    "only one external source"), "`robust_sd` is ignored because robust_weight = 0")
})
