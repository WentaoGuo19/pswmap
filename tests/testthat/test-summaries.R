test_that("Hajek mean, design-based variance and Rubin's rules are computed as defined", {
  set.seed(2)
  w = psw_weights(ex_trial(), ex_ext(), mi_m = 4)
  s = psw_summaries(w)
  src = w$sources$H1
  per = apply(src$w, 2, function(wi) { m = sum(wi * src$y) / sum(wi); c(m, sum(wi^2 * (src$y - m)^2) / sum(wi)^2) })
  expect_equal(s$mean[1], mean(per[1, ]))
  expect_equal(s$se[1], sqrt(mean(per[2, ]) + (1 + 1 / 4) * stats::var(per[1, ])))
  expect_equal(s$kish_ess[1], mean(apply(src$w, 2, function(wi) sum(wi)^2 / sum(wi^2))))
})
