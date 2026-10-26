#' Example data: a randomized trial and three external control sources
#'
#' Simulated data for trying out the package. A randomized trial is accompanied by three external control sources
#' whose covariate distributions differ from the trial's. Five covariates are partly missing in the external data,
#' while the trial is complete. The data follow the design of the simulation study of Guo and Zhu:
#' - low covariate shift;
#' - no unmeasured differences in the outcome between the sources and the trial;
#' - a true treatment effect of 0.7;
#' - covariates `X2`, `X4`, `X6`, `X8` and `X10` missing at random: the probability that `X2` is missing depends on
#'   the always-observed `X1`, that of `X4` on `X3`, and so on. About 20% of external subjects have at least one
#'   missing covariate.
#'
#' The outcome is generated as \eqn{Y = \beta_0^\top X + 0.7\,T + \varepsilon}{Y = beta_0' X + 0.7 T + epsilon} in the trial and
#' \eqn{Y = \beta_0^\top X + \varepsilon}{Y = beta_0' X + epsilon} in the external sources, where \eqn{\varepsilon \sim N(0, 1)}{epsilon ~ N(0, 1)}. The
#' coefficients \eqn{\beta_0}{beta_0} are 0.1 for the always-observed covariates and 0.3 for the covariates that may be
#' missing.
#'
#' @format A list with two elements:
#' \describe{
#'   \item{`trial`}{a data frame of 150 subjects: the outcome `Y`, the treatment `trt` (100 treated, 50 control) and
#'     the covariates `X1`-`X10` (`X1`-`X4` binary, `X5`-`X10` continuous).}
#'   \item{`external`}{a named list of three data frames, `H1`, `H2` and `H3`, of 200 external controls each, with
#'     `Y` and `X1`-`X10`.}
#' }
#' @source Simulated with the data-generating functions of the simulation study (`R/simulate_data.R` of the
#'   pswmap-paper repository, <https://github.com/WentaoGuo19/pswmap-paper>). The script that creates the data is
#'   `data-raw/pswmap_example.R` in the source repository of pswmap.
#' @examples
#' data(pswmap_example)
#' str(pswmap_example$trial)
#' # share of external subjects with at least one missing covariate
#' sapply(pswmap_example$external, function(d) mean(!complete.cases(d)))
"pswmap_example"
