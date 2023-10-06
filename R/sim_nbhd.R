#' Simulate a neighborhood for a respondent
#'
#' @param fit the model fit, from [neighborhood_model()]
#' @param new_resp a single-row respondent data frame to make predictions from.
#'  Should have a `neighborhood` column with codes matching `block_gr$blocks`.
#' @param draws the number of simulated neighborhoods per respondent
#' @param resp_id the resondent ID, used to get the random effects. If `NULL`
#'   a new random effect will be simulated for each draw. If `NA`, the random
#'   effect will be set to zero. If a numeric, the random effect will be set to
#'   this value.
#' @param block_d the block data frame. Should have a `centroid` geography
#'   column.
#' @param block_gr An adjacency graph object: a list containing an element
#'   `graph` with the adjacency list, and `blocks` a character vector of GEOIDs
#'   or codes corresponding to the indices in `graph`.
#' @param max_ring the maximum graph distance from the starting block to allow.
#' @param proc_fn a processing function that is used to prepare raw model data
#'   for fitting
#'
#' @returns A list of `draws` integer vectors containing the indices of the
#'   blocks (in `block_d`) making up each simulated neighborhood.
#'
#' @export
simulate_neighborhood = function(fit, new_resp, draws=1, resp_id=NULL,
                                 block_d, block_gr, max_ring=10L, proc_fn=function(x) x) {
    stopifnot(inherits(fit, "nbhd_fit"))
    draw_idxs = sample.int(fit$draws, draws, replace=TRUE)
    coefs = t(cbind(draws_of(fit$post$coefs * fit$post$alpha),
                    draws_of(fit$post$alpha))[draw_idxs, , drop=FALSE])
    # reorder to match model matrix
    coefs = coefs[match(colnames(fit$X), c(names(fit$post$coef), "log(dist)")),
                  , drop=FALSE]
    # figure out random effects
    if (is.null(resp_id)) {
        ranef = rnorm(draws, 0, draws_of(fit$post$id_sd)[draw_idxs, ])
        resp_id = -1L
    } else if (is.numeric(resp_id) && !is.integer(resp_id)) {
        ranef = rep(resp_id, draws)
        resp_id = -1L
    } else if (as.character(resp_id) %in% names(fit$ids)) {
        id_idx = fit$ids[as.character(resp_id) == names(fit$ids)]
        ranef = draws_of(fit$post$id_ranef[id_idx])[draw_idxs, ]
    } else if (is.na(resp_id)) {
        ranef = rep(0.0, draws)
        resp_id = -1L
    } else {
        warning("`resp_id` not found; setting random effect to zero.")
        ranef = rep(0.0, draws)
        resp_id = -1L
    }

    start_fips = new_resp$neighborhood[[1]][1]
    start_idx = match(start_fips, block_gr$blocks)
    stopifnot(!is.na(start_idx))
    new_resp = select(new_resp, -"neighborhood")
    new_resp$id = as.character(resp_id)

    centr_col = which(names(block_d) == "centroid")
    if (!inherits(block_d$centroid, "s2_geography"))
        block_d$centroid = as_s2_geography(block_d$centroid)

    done = FALSE
    incl = matrix(rep(TRUE, draws), nrow=1, ncol=draws) # the output: neighborhood inclusion
    last_ring = 1L # how many blocks in last ring
    last_dists = 0
    rownames(incl) = start_idx
    ring_range = c(0L, max_ring)
    active_cols = rep(TRUE, draws) # which neighborhoods are still being drawn
    while (!done) {
        relevant = get_within_ring(ring_range[2], start_idx, block_gr$graph)
        idxs = relevant$idx[relevant$rings > ring_range[1]]
        rings = relevant$rings[relevant$rings > ring_range[1]]
        dists = as.numeric(s2_distance_matrix(block_d$centroid[start_idx],
                                              block_d$centroid[idxs]))
        order_idx = order(rings, dists)
        idxs = idxs[order_idx]

        pred_d = bind_cols(block_d[idxs, -centr_col], new_resp)
        pred_d$dist = dists[order_idx]
        pred_d$ring = rings[order_idx]
        pred_d$incl = 0L
        pred_d = proc_fn(pred_d)

        X = model.matrix(fit$terms, data=pred_d, xlev=fit$xlevels)
        if (nrow(X) != length(idxs)) stop("Missing data in prediciton frame.")
        new_incl = sim_incl(X %*% coefs[, active_cols, drop=FALSE],
                            ranef[active_cols])
        old_idxs = seq(nrow(incl)-last_ring+1L, nrow(incl))
        new_incl = rbind(incl[old_idxs, active_cols, drop=F], new_incl)

        check_idxs = c(as.integer(rownames(incl)[old_idxs]), idxs)
        check_ring = c(rep(ring_range[1], last_ring), pred_d$ring)
        check_dists = c(last_dists, pred_d$dist)
        new_incl = fix_incl_mat(check_idxs, new_incl,
                                check_ring, check_dists, block_gr$graph)
        rownames(new_incl) = check_idxs

        add_incl = matrix(FALSE, nrow=nrow(new_incl)-last_ring, ncol=draws)
        add_incl[, active_cols] = new_incl[-1:-last_ring, ]
        rownames(add_incl) = rownames(new_incl)[-1:-last_ring]
        incl = rbind(incl, add_incl)
        last_ring = sum(pred_d$ring == ring_range[2])
        last_dists = pred_d$dist[pred_d$ring == ring_range[2]]

        active_cols[active_cols] = apply(new_incl[pred_d$ring == ring_range[2],
                                                  , drop=FALSE], 2, any)
        if (!any(active_cols)) {
            done = TRUE
        } else {
            ring_range = ring_range + diff(ring_range)
        }
    }

    blocks = as.integer(rownames(incl))
    lapply(1:draws, function(i) blocks[incl[, i]])
}
