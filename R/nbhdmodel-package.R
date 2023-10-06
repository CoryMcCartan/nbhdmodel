#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang env enexpr expr caller_env
#' @import dplyr
#' @import ggplot2
#' @import sf
#' @import s2
#' @import rstan
#' @importFrom posterior quantile2 draws_of
#' @importFrom stringr str_detect str_replace str_sub str_length str_c
#' @importFrom fastmatch fmatch fmatch.hash
#' @importFrom Rcpp evalCpp
#' @useDynLib nbhdmodel, .registration = TRUE
## usethis namespace: end
NULL

sf::sf_use_s2(TRUE)
