#' Get Posterior Block Inclusion Probabilities
#'
#' @param fit the model fit, from [neighborhood_model()]
#' @param new_resp a single-row respondent data frame to make predictions from
#' @param resp_id the resondent ID, used to get the random effects. If `NULL`
#'   a new random effect will be simulated for each draw. If `NA`, the random
#'   effect will be set to zero.
#' @param block_d the block data frame
#' @param proc_fn a processing function that is used to prepare raw model data
#'   for fitting
#' @param use_distance if `FALSE`, remove block-to-home distance from the linear
#'   predictor.
#'
#' @returns a numeric vector of inclusion probabilities
#'
#' @export
post_incl = function(fit, new_resp, resp_id=NULL, block_d,
                     proc_fn=function(x) x, use_distance=TRUE) {
    stopifnot(inherits(fit, "nbhd_fit"))
    coefs = t(cbind(draws_of(fit$post$coefs * fit$post$alpha),
                    draws_of(fit$post$alpha)))
    # reorder to match model matrix
    coefs = coefs[match(colnames(fit$X), c(names(fit$post$coef), "log(dist)")),
                  , drop=FALSE]
    # figure out random effects
    draws = fit$draws
    if (is.null(resp_id)) {
        ranef = rnorm(draws, 0, draws_of(fit$post$id_sd))
        resp_id = -1L
    } else if (is.numeric(resp_id) && !is.integer(resp_id)) {
        ranef = rep(resp_id, draws)
        resp_id = -1L
    } else if (as.character(resp_id) %in% names(fit$ids)) {
        id_idx = fit$ids[as.character(resp_id) == names(fit$ids)]
        ranef = draws_of(fit$post$id_ranef[id_idx])
    } else if (is.na(resp_id)) {
        ranef = rep(0.0, draws)
        resp_id = -1L
    } else {
        warning("`resp_id` not found; setting random effect to zero.")
        ranef = rep(0.0, draws)
        resp_id = -1L
    }

    centr_col = which(names(block_d) == "centroid")
    start_fips = new_resp$neighborhood[[1]][1]
    start_idx = match(start_fips, block_d$fips)
    pred_d = bind_cols(block_d[, -centr_col], new_resp)
    if (use_distance) {
        dists = as.numeric(s2_distance_matrix(block_d$centroid[start_idx],
                                              block_d$centroid))
        pred_d$dist = dists
    } else {
        pred_d$dist = 0
    }
    pred_d$incl = 0
    pred_d = proc_fn(pred_d)

    X = model.matrix(fit$terms, data=pred_d, xlev=fit$xlevels)
    if (use_distance) {
        pr_incl(X %*% coefs, ranef)
    } else {
        dist_col = which(colnames(X) == "log(dist)")
        pr_incl(X[, -dist_col] %*% coefs[-dist_col, ], ranef)
    }
}

#' Get Posterior Mean of Effective Block Distance
#'
#' Not exported. Assumes `block_d` has column `fips` which is used inside
#' the `neighborhood` column of `new_resp`..
#'
#' @param fit the model fit, from [neighborhood_model()]
#' @param new_resp a single-row respondent data frame to make predictions from
#' @param block_d the block data frame
#' @param proc_fn a processing function that is used to prepare raw model data
#'   for fitting
#'
#' @returns a numeric vector of effective distances
#'
#' @export
eff_dist = function(fit, new_resp, block_d, proc_fn=function(x) x) {
    stopifnot(inherits(fit, "nbhd_fit"))
    coefs = t(cbind(draws_of(fit$post$coefs * fit$post$alpha),
                    draws_of(fit$post$alpha)))
    # reorder to match model matrix
    coefs = coefs[match(colnames(fit$X), c(names(fit$post$coef), "log(dist)")),
                  , drop=FALSE]

    centr_col = which(names(block_d) == "centroid")
    start_fips = new_resp$neighborhood[[1]][1]
    start_idx = match(start_fips, block_d$fips)
    pred_d = bind_cols(block_d[, -centr_col], new_resp)
    pred_d$dist = 0
    pred_d$incl = 0
    pred_d = proc_fn(pred_d)

    X = model.matrix(fit$terms, data=pred_d, xlev=fit$xlevels)
    drop_coef = which(colnames(X) %in% c("(Intercept)", "log(dist)"))
    rowMeans(exp(X[, -drop_coef] %*% coefs[-drop_coef, ]))
}

