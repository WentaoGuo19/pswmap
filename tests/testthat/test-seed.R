test_that("the same seed gives the same analysis, and a different seed a different one", {
  run = function(seed) pswmap(ex_trial(), ex_ext(), tau_c = 16, mi_m = 2, map = fast_map(), seed = seed)
  a = run(11); b = run(11); c = run(12)
  expect_identical(a$effect, b$effect)
  expect_identical(a$sources, b$sources)
  expect_false(identical(a$effect, c$effect))
})

test_that("seed leaves the caller's random number stream unchanged", {
  set.seed(99); x_ref = runif(3)
  set.seed(99)
  pswmap(ex_trial(), borrow = FALSE, seed = 1)
  expect_identical(runif(3), x_ref)
})

test_that("seed is equivalent to set.seed() before the call", {
  a = pswmap(ex_trial(), ex_ext(), tau_c = 16, missing = "cc", map = fast_map(), seed = 5)
  set.seed(5)
  b = pswmap(ex_trial(), ex_ext(), tau_c = 16, missing = "cc", map = fast_map())
  expect_identical(a$effect, b$effect)
})

test_that("seed must be a single number", {
  expect_error(pswmap(ex_trial(), borrow = FALSE, seed = "a"), "seed")
  expect_error(pswmap(ex_trial(), borrow = FALSE, seed = 1:2), "seed")
})
