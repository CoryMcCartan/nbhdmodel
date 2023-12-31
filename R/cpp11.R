# Generated by cpp11: do not edit by hand

find_relevant <- function(nbhd, graph) {
  .Call(`_nbhdmodel_find_relevant`, nbhd, graph)
}

check_incl <- function(idx, incl, rings, dists, graph) {
  .Call(`_nbhdmodel_check_incl`, idx, incl, rings, dists, graph)
}

pr_incl <- function(base_linpred, ranef) {
  .Call(`_nbhdmodel_pr_incl`, base_linpred, ranef)
}

get_within_ring <- function(ring, start, graph) {
  .Call(`_nbhdmodel_get_within_ring`, ring, start, graph)
}

sim_incl <- function(base_linpred, ranef) {
  .Call(`_nbhdmodel_sim_incl`, base_linpred, ranef)
}

fix_incl_mat <- function(idx, incl, rings, dists, graph) {
  .Call(`_nbhdmodel_fix_incl_mat`, idx, incl, rings, dists, graph)
}
