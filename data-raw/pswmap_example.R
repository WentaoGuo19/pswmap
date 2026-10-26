# Generates data/pswmap_example.rda with the data-generating functions of the simulation study of Guo and Zhu
# (R/simulate_data.R of the pswmap-paper repository, https://github.com/WentaoGuo19/pswmap-paper): low covariate
# shift, linear selection, no outcome heterogeneity, treatment effect 0.7, and covariates X2, X4, X6, X8 and X10
# missing at random in the external data (20% of external subjects with at least one missing covariate).
#
# Run from the package root. Set `paper_repo` to your copy of the pswmap-paper repository (it needs the package MASS).
paper_repo = "../pswmap-paper"
source(file.path(paper_repo, "R", "simulate_data.R"))
set.seed(2026)
dc = simdat_complete(n_treat = 100, n_control = 50, n_hist = rep(200, 3), mu_current = 1, mu_hist = c(1, 0.5, 1.2),
                     delta = c(0, 0, 0), Vtimer = c(1, 1, 1), beta1 = 0.7, beta0_com = 0.1, beta0_mis = 0.3)
dm = simdat_missing(dc, mechanism = "MAR", row_rate = 0.2)
xv = paste0("X", 1:10)
trial = data.frame(Y = dm$RCT$Y, trt = dm$RCT$treatment, dm$RCT[, xv])
external = lapply(dm$RWD, function(d) data.frame(Y = d$Y, d[, xv]))
names(external) = paste0("H", seq_along(external))
pswmap_example = list(trial = trial, external = external)
dir.create("data", showWarnings = FALSE)
save(pswmap_example, file = "data/pswmap_example.rda", compress = "xz")
