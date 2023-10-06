
#' Create a model data frame for an individual respondent
#'
#' This function combines individual, neighborhood, and geographic information
#' to produce a data frame suitable for use in fitting a [neighborhood_model()].
#' Can be applied in a loop over respondents to generate a model frame for
#' an entire sample.
#'
#' @param row A single data frame row containing relevant individual covariates
#'   for a respondent.
#' @param nbhd The respondent's neighborhood as a vector of indices indexing the
#'   blocks in `block_gr`.
#' @param block_d A data frame of census blocks, including a column `centroid`
#'   with an `s2` point geography of each block's centroid.
#' @param block_gr An adjacency graph object: a list containing an element
#'   `graph` with the adjacency list, and `blocks` a character vector of GEOIDs
#'   or codes corresponding to the indices in `graph`.
#'
#' @returns A tibble that can be used inside a modeling function. Will contain
#' the entries in `row`, plus the relevant entries in `block_d` for each
#' block, plus columns:
#'
#' - `ring` containing the "ring" indicator around the residence: 0 indicates
#' the respondent's block, 1 for blocks touching the residential block, 2 for blocks
#' touching those, etc.
#' - `incl` a binary indicator for whether the block is in the neighborhood
#' - `dist` the distance to the respondent's block
#' - `frac_con` the fraction of nearer blocks in the neighborhood this block
#' is connected to.
#'
#' @export
calc_indiv_frame = function(row, nbhd, block_d, block_gr) {
    if (length(nbhd) <= 1) stop("Neighborhood must have at least two blocks.")

    rel = find_relevant(nbhd, block_gr$graph)
    dists = as.numeric(s2_distance_matrix(block_d$centroid[nbhd[1]],
                                          block_d$centroid[rel$idx]))
    frac_con = check_incl(rel$idx, rel$incl, rel$rings, dists, block_gr$graph)
    frac_con[1] = 1
    keep_idx = frac_con > 0

    block_d$centroid = NULL
    out = bind_cols(block_d[rel$idx[keep_idx],], row)
    out$ring = rel$rings[keep_idx]
    out$incl = rel$incl[keep_idx]
    out$dist = dists[keep_idx]
    out$frac_con = frac_con[keep_idx]
    order_idx = order(out$ring, out$dist)
    out[order_idx,]
}
