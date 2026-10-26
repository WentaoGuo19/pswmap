# small example data and a fast MCMC setting for the tests
ex_trial = function() pswmap_example$trial
ex_ext = function() pswmap_example$external
fast_map = function() map_control(iter = 4000, warmup = 1000, chains = 2, thin = 1)
