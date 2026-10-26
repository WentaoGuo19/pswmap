# Propensity score learners, missing-data strategies and imputation (internal)

# Default tuning of each learner (the settings recommended by Guo and Zhu); the user's learner_args override them.
learner_defaults = function(learner) {
  switch(learner,
         logit = list(),
         bart = list(),   #bartMachine's own defaults
         rf = list(num.trees = 500, num.threads = 1),
         xgb = list(eta = 0.1, max_depth = 3, min_child_weight = 1, nrounds = 100, nfolds = 5, nthread = 1))
}

need_pkg = function(pkg, learner) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("learner = \"", learner, "\" needs the package ", pkg, ": install.packages(\"", pkg, "\")", call. = FALSE)
}


# Fit one propensity score learner: P(Z = 1 | X) for the stacked data X with trial indicator z (1 trial, 0 external).
# native = TRUE lets the learner use the missing values (missingness indicators for logit, BARTm for bart, ranger's
# na.learn for rf, XGBoost's default split direction for xgb); native = FALSE requires complete X.
# Returns fitted probabilities: in-sample (logit, bart), out-of-bag (rf) or cross-fitted (xgb).
fit_learner = function(X, z, learner, native = FALSE, args = list()) {
  if (!native && anyNA(X)) stop("missing covariates reached a learner without missing-data handling.", call. = FALSE)
  args = utils::modifyList(learner_defaults(learner), args)

  if (learner == "logit") {
    if (native && anyNA(X)) {
      for (v in names(X)[colSums(is.na(X)) > 0]) {
        X[[paste0("M_", v)]] = as.numeric(is.na(X[[v]]))
        X[[v]][is.na(X[[v]])] = mean(X[[v]], na.rm = TRUE)
      }
    }
    #missingness indicators separate incomplete external subjects from the (complete) trial: fitted probabilities
    #near 0 and a separation warning are expected
    fit = suppressWarnings(stats::glm(z ~ ., family = stats::binomial, data = cbind(X, z = z)))
    return(unname(stats::fitted(fit)))
  }

  if (learner == "bart") {
    need_pkg("bartMachine", learner)
    fit = do.call(bartMachine::bartMachine,
                  c(list(X = X, y = factor(z, levels = c(1, 0)), use_missing_data = native,
                         use_missing_data_dummies_as_covars = native, verbose = FALSE), args))
    return(fit$p_hat_train)   #posterior mean probability of the first level (trial)
  }

  if (learner == "rf") {
    need_pkg("ranger", learner)
    fit = do.call(ranger::ranger, c(list(x = X, y = factor(z, levels = c(0, 1)), probability = TRUE,
                                         na.action = "na.learn"), args))
    p = fit$predictions[, "1"]
    #subjects never out-of-bag (rare): in-sample prediction
    if (anyNA(p)) p[is.na(p)] = stats::predict(fit, X[is.na(p), , drop = FALSE])$predictions[, "1"]
    return(p)
  }

  if (learner == "xgb") {
    need_pkg("xgboost", learner)
    nfolds = args$nfolds; nrounds = args$nrounds
    params = c(list(objective = "binary:logistic"), args[setdiff(names(args), c("nfolds", "nrounds"))])
    folds = integer(length(z))
    for (g in 0:1) folds[z == g] = sample(rep_len(seq_len(nfolds), sum(z == g)))   #stratified by trial indicator
    Xm = as.matrix(X)
    p = numeric(length(z))
    for (f in seq_len(nfolds)) {
      tr = folds != f
      dtr = xgboost::xgb.DMatrix(Xm[tr, , drop = FALSE], label = z[tr], missing = NA)
      fit = xgboost::xgb.train(params = params, data = dtr, nrounds = nrounds, verbose = 0)
      p[!tr] = stats::predict(fit, xgboost::xgb.DMatrix(Xm[!tr, , drop = FALSE], missing = NA))
    }
    return(p)
  }
  stop("`learner` must be one of \"logit\", \"bart\", \"rf\", \"xgb\".", call. = FALSE)
}


# Propensity scores of the external subjects for one completed (or incomplete, for native) set of covariates.
# trial_X: complete trial covariates; ext_X: list of external covariate data frames.
# strategy: "full" (complete data), "cc" (complete-case: incomplete subjects get PS 0) or "native".
# scope: "source" (one model per source: trial vs that source) or "pooled" (one model: trial vs all sources).
# Returns a list (one numeric vector per source) of PS for the external subjects.
ps_external = function(trial_X, ext_X, strategy, scope, learner, args) {
  fit_one = function(Xe_list) {
    #fit trial vs the stacked sources in Xe_list; return PS of the external rows, split by source
    Xe = do.call(rbind, Xe_list)
    n_e = vapply(Xe_list, nrow, integer(1))
    ps = numeric(nrow(Xe))
    if (strategy == "cc") {
      cc = stats::complete.cases(Xe)
      X = rbind(trial_X, Xe[cc, , drop = FALSE])
      z = c(rep(1, nrow(trial_X)), rep(0, sum(cc)))
      ps[cc] = fit_learner(X, z, learner, native = FALSE, args = args)[z == 0]
    } else {
      X = rbind(trial_X, Xe)
      z = c(rep(1, nrow(trial_X)), rep(0, nrow(Xe)))
      ps = fit_learner(X, z, learner, native = strategy == "native", args = args)[z == 0]
    }
    split(ps, rep(seq_along(Xe_list), n_e))
  }
  if (scope == "source") {
    out = lapply(ext_X, function(Xk) fit_one(list(Xk))[[1]])
  } else {
    out = fit_one(ext_X)
  }
  names(out) = names(ext_X)
  out
}


# Multiply impute the missing covariates of one external source with mice, using that source's subjects only and
# with the outcome in the imputation model. engine "default": logistic regression for binary 0/1 covariates and
# predictive mean matching for the others; "cart" or "rf": that method for every incomplete covariate.
# Returns a list of m completed covariate data frames.
impute_source = function(d, covariates, outcome, m, engine, mi_args) {
  if (!anyNA(d[, covariates])) return(rep(list(d[, covariates, drop = FALSE]), m))
  d = d[, c(covariates, outcome), drop = FALSE]
  meth = mice::make.method(d)
  inc = names(meth)[meth != ""]
  if (engine == "default") {
    bin = inc[vapply(inc, function(v) all(stats::na.omit(d[[v]]) %in% c(0, 1)), logical(1))]
    meth[bin] = "logreg"
  } else if (engine %in% c("cart", "rf")) {
    bin = character(0)
    meth[inc] = engine
  } else stop("`mi_engine` must be one of \"default\", \"cart\", \"rf\".", call. = FALSE)
  #logistic-regression imputation needs factors: code binary covariates as factors, back to 0/1 afterwards
  for (v in bin) d[[v]] = factor(d[[v]], levels = c(0, 1))
  mids = suppressWarnings(do.call(mice::mice, c(list(data = d, m = m, method = meth, printFlag = FALSE), mi_args)))
  lapply(seq_len(m), function(i) {
    ci = mice::complete(mids, action = i)
    for (v in bin) ci[[v]] = as.numeric(as.character(ci[[v]]))
    if (anyNA(ci[, covariates]))
      stop("mice left missing values in external covariates; try another `mi_engine` or `mi_args`.", call. = FALSE)
    ci[, covariates, drop = FALSE]
  })
}
