test_that("complete-case weighting gives incomplete subjects weight zero", {
  w = psw_weights(ex_trial(), ex_ext(), missing = "cc")
  expect_equal(w$missing, "cc")
  for (s in w$sources) expect_true(all(s$w[s$incomplete, ] == 0))
  expect_equal(psw_summaries(w)$relweight_incomplete, c(0, 0, 0))
})

test_that("logistic regression with missingness indicators converges to complete-case weighting", {
  w = psw_weights(ex_trial(), ex_ext(), missing = "native")
  expect_true(all(psw_summaries(w)$relweight_incomplete < 1e-3))
})

test_that("multiple imputation keeps incomplete subjects and gives one weight column per imputation", {
  set.seed(1)
  w = psw_weights(ex_trial(), ex_ext(), mi_m = 3)
  expect_equal(w$missing, "mi")
  expect_equal(w$m, 3)
  for (s in w$sources) expect_equal(ncol(s$w), 3)
  expect_true(all(psw_summaries(w)$relweight_incomplete > 0.5))
})

test_that("complete external data use no missing-data strategy", {
  ext = lapply(ex_ext(), function(d) d[stats::complete.cases(d), ])
  w = psw_weights(ex_trial(), ext, missing = "cc")
  expect_equal(w$missing, "none")
  expect_equal(w$m, 1)
  expect_true(all(is.na(psw_summaries(w)$relweight_incomplete)))
})

test_that("the propensity score cap bounds the weights", {
  w = psw_weights(ex_trial(), ex_ext(), missing = "cc", ps_cap = 0.9)
  expect_true(max(unlist(lapply(w$sources, `[[`, "w"))) <= 0.9 / 0.1 + 1e-12)
  expect_error(psw_weights(ex_trial(), ex_ext(), ps_cap = 1.2), "ps_cap")
})

test_that("a pooled propensity score triggers a warning", {
  expect_warning(psw_weights(ex_trial(), ex_ext(), missing = "cc", ps_scope = "pooled"), "not recommended")
})

test_that("the flexible learners run", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("xgboost")
  set.seed(1)
  for (l in c("rf", "xgb")) {
    w = psw_weights(ex_trial(), ex_ext(), missing = "native", learner = l)
    r = psw_summaries(w)$relweight_incomplete
    expect_true(all(r > 0 & r < 1))
  }
})
