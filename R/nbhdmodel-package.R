#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang env enexpr expr caller_env
#' @import dplyr
#' @import sf
#' @import s2
#' @importFrom posterior quantile2 draws_of
#' @importFrom stringr str_detect str_replace str_sub str_length str_c str_starts str_split_fixed
#' @importFrom fastmatch fmatch fmatch.hash
#' @importFrom stats .getXlevels fitted median model.frame model.matrix
#' @importFrom stats resid rnorm terms update.formula
#' @importFrom utils capture.output globalVariables
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel CxxFlags
#' @useDynLib nbhdmodel, .registration = TRUE
## usethis namespace: end
NULL

sf::sf_use_s2(TRUE)
globalVariables(c("nm", ".", "ring", "variable", "med", "out_l", "out_h",
                  "in_l", "in_h", "rhat", "ess_tail"))
