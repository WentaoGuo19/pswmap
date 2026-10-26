test_that("pswmap borrows and returns a robust prior and a posterior", {
  set.seed(3)
  fit = pswmap(ex_trial(), ex_ext(), tau_c = 16, missing = "cc", map = fast_map())
  expect_s3_class(fit, "pswmap")
  expect_true(fit$borrow)
  expect_equal(nrow(fit$sources), 3)
  expect_true(fit$prior$ess > 0)
  expect_equal(unname(fit$prior$robust[1, ncol(fit$prior$robust)]), 0.2)   # robust weight
  expect_true(fit$effect["lower"] < fit$effect["upper"])
  expect_true(fit$prob_positive >= 0 && fit$prob_positive <= 1)
  expect_output(print(fit), "PSW-MAP analysis")
  expect_output(print(summary(fit)), "Robust PSW-MAP prior")
})

test_that("borrow = FALSE analyses the trial alone, with or without external data", {
  set.seed(4)
  a = pswmap(ex_trial(), borrow = FALSE)
  expect_false(a$borrow)
  expect_null(a$prior)
  tr = ex_trial()
  expect_equal(unname(a$effect["mean"]), mean(tr$Y[tr$trt == 1]) - mean(tr$Y[tr$trt == 0]), tolerance = 0.02)
  expect_output(print(a), "trial alone")
})

test_that("with the same seed, different values of tau_c reuse the same weights", {
  a = pswmap(ex_trial(), ex_ext(), tau_c = 32, missing = "cc", map = fast_map(), seed = 5)
  b = pswmap(ex_trial(), ex_ext(), tau_c = 4, missing = "cc", map = fast_map(), seed = 5)
  expect_identical(a$sources, b$sources)
  expect_true(a$prior$ess > b$prior$ess)   # a smaller c borrows less
})

test_that("robust_weight = 0 returns the MAP prior unchanged", {
  set.seed(6)
  s = psw_summaries(psw_weights(ex_trial(), ex_ext(), missing = "cc"))
  p = psw_map_prior(s, sigma_ref = 1.3, tau_c = 16, control = fast_map(), robust_weight = 0)
  expect_identical(p$robust, p$map)
  expect_error(psw_map_prior(s, sigma_ref = 1.3), "tau_c")
})
