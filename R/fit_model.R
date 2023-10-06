#' Fit the Neighborhood Model
#'
#' Fits the neighborhood GLMM using a provided formula and using a
#' Bernoulli outcome with cloglog link.
#'
#' @param formula a one-sided model formula. The actual formula used to fit the
#'   GLMM will be generated from this one: it will contain a term for log(dist)
#'   (where `dist` must be the column that encodes distance), and will have random
#'   effects based on the `id` column.  The column `incl` will be used as the
#'   left-hand-side variable indicating that a block is in the neighborhood.
#' @param data the model data frame.
#'   Should have an `incl` column for block inclusion, an  `id` column with
#'   respondent IDs, and a `dist` column with distances.
#' @param prior_coef_scale the scale of the prior on the standardized predictors.
#' @param draws the number of approximate posterior draws to generate
#' @param imp_samp whether to perform importance resampling on the approximate draws.
#' @param init initial values for the model fitting function [rstan::optimizing()].
#' @param ... other arguments to to the model fitting function [rstan::optimizing()].
#' @param hessian whether to compute the Hessian. Required for full inference.
#' @param verbose if `TRUE`, show verbose optimization output.
#'
#' @returns a fitted model object of class [nbhd_fit], which is a list which
#' includes the following elements:
#'
#' - `map`, the MAP estimates for the parameters.
#' - `vcov`, the covariance matrix for the MAP estimates, calculated from the
#' Hessian of the log posterior.
#' - `raw_ids` the vector of ids
#' - `X` the design matrix
#' - `y` the outcome value (1 - `incl`)
#' - `post`, the approximate posterior samples, from [posterior::draws_rvars]
#' - `lp`, the log posterior probability of each sample
#' - `lp_norm`, the log probability of Normal approx. to posterior for each sample
#'
#' @export
neighborhood_model = function(formula, data, prior_coef_scale=2.5, draws=1000,
                              imp_samp=TRUE, init=0, ..., hessian=TRUE, verbose=FALSE) {
    if (any(is.na(data)))
        stop("Data must not contain NA values.")
    fit_d = filter(data, ring > 0) %>%
        mutate(across(where(is.factor), droplevels))
    fit_form = update.formula(formula, ~ . + log(dist))

    X_fr = model.frame(fit_form, data=fit_d, drop.unused.levels=FALSE)
    X = model.matrix(fit_form, data=fit_d, drop.unused.levels=FALSE)
    id = as.integer(as.factor(fit_d$id))
    stan_d = list(
        N = nrow(fit_d),
        Y = 1L - fit_d$incl,
        X = X,
        K = ncol(X),
        id = id,
        N_id = length(unique(id)),
        bQ_prior_scale = prior_coef_scale,
        bQ_prior_df = 20L,
        prior_only = 0
    )

    if (qr(stan_d$X)$rank < ncol(stan_d$X)) {
        stop("Covariate matrix is not full rank.")
    }

    fit = rstan::optimizing(stanmodels$nbhd_glmm,
                            data=stan_d, init=init, hessian=hessian, draws=draws,
                            verbose=verbose, importance_resampling=imp_samp, ...)
    if (!hessian) return(fit)
    if (fit$return_code != 0) cli::cli_abort("Finding MAP failed.")

    drop_pars = which(str_detect(names(fit$par), "bQ\\[") |
                      names(fit$par) == "Intercept")

    #vcov = - chol2inv(R[-drop_pars, -drop_pars])
    vcov = -solve(fit$hessian[-drop_pars, -drop_pars])

    new_names = names(fit$par)
    alpha_idx = which(colnames(stan_d$X) == "log(dist)")-1L
    new_names = str_replace(new_names, "z_1\\[", "id_ranef[")
    new_names = str_replace(new_names, "sd_1", "id_sd")
    new_names = str_replace(new_names, paste0("b\\[", alpha_idx, "]"), "alpha")
    new_names = str_replace(new_names, "b_Intercept", "b[0]")
    new_names = str_replace(new_names, "b\\[", "coefs[")
    names(fit$par) = new_names
    colnames(fit$theta_tilde) = new_names

    n_pars = length(new_names[-drop_pars])
    fixef_start = which.max(str_detect(new_names[-drop_pars], "coefs\\[1"))
    fixef_idx = c(n_pars, fixef_start:(n_pars-1L))
    ranef_idx = 2:(fixef_start-1L)
    names(fixef_idx) = colnames(stan_d$X)
    names(ranef_idx) = levels(as.factor(fit_d$id))

    ids_out = seq_len(stan_d$N_id)
    names(ids_out) = levels(as.factor(fit_d$id))
    out = structure(
        list(map = fit$par[-drop_pars],
             vcov = vcov,
             post = posterior::as_draws_rvars(fit$theta_tilde[, -drop_pars]),
             fixef = fixef_idx,
             ranef = ranef_idx,
             lp = fit$log_p,
             lp_norm = fit$log_g,
             raw_ids = stan_d$id,
             ids = ids_out,
             y = stan_d$Y,
             X = stan_d$X,
             draws = draws,
             call = sys.call(),
             terms = terms(fit_form),
             xlevels = .getXlevels(terms(fit_form), X_fr),
             prior = list(coef_scale = prior_coef_scale)
             ),
        class = "nbhd_fit",
        imp_samp = imp_samp)
    names(out$post$id_ranef) = levels(as.factor(fit_d$id))
    names(out$post$coefs) = colnames(stan_d$X)[-(alpha_idx+1L)]
    out$post$coefs = out$post$coefs / out$post$alpha
    out$post$id_ranef = out$post$id_ranef * out$post$id_sd

    out
}

#' @export
print.nbhd_fit = function(x, ...) {
    cli::cat_line("Neighborhood Generalized Linear Mixed Model\n\nCall:")
    cli::cat_line(rlang::expr_text(x$call))
    cli::cat_line("\nRespondents: ", length(x$ids))
    cli::cat_line("Blocks:      ", nrow(x$X))
    cli::cat_line("\nIndividual intercept std. dev.: ",
                  str_sub(capture.output(print(x$post$id_sd))[2], 5))
    cli::cat_line("Kernel shape parameter: ",
                  str_sub(capture.output(print(x$post$alpha))[2], 5))
    cli::cat_line("\nCoefficients:")
    cli::cat_line(as.list(capture.output(print(x$post$coefs))[-1]))
}

#' Functions for working with neighborhood fits
#'
#' @name nbhd_fit
NULL

#' @param object,x a `nbhd_fit` object
#' @param ... Ignored.
#'
#' @method summary nbhd_fit
#' @rdname nbhd_fit
#' @export
summary.nbhd_fit = function(object, ...) {
    summary(object$post[c("coefs", "alpha")]) %>%
        mutate(variable = if_else(str_starts(variable, "coefs"),
                                  str_sub(variable, 7, -2),
                                  variable)) %>%
        select(-"rhat":-"ess_tail")
}

#' @method coef nbhd_fit
#' @rdname nbhd_fit
#' @export
coef.nbhd_fit = function(object, ...) {
    coefs = object$map[object$fixef]
    names(coefs) = names(object$fixef)
    coefs
}


#' @importFrom lme4 fixef ranef
#' @export
lme4::fixef
#' @export
lme4::ranef

#' @method fixef nbhd_fit
#' @rdname nbhd_fit
#' @export
fixef.nbhd_fit = function(object, ...) {
    coef.nbhd_fit(object)
}
#' @method ranef nbhd_fit
#' @rdname nbhd_fit
#' @export
ranef.nbhd_fit = function(object, ...) {
    ranefs = object$map[object$ranef]
    names(ranefs) = names(object$ranef)
    ranefs
}

#' @method fitted nbhd_fit
#' @rdname nbhd_fit
#' @export
fitted.nbhd_fit = function(object, ...) {
    linpred = object$X %*% coef.nbhd_fit(object) + object$map[1]*ranef(object)[object$raw_ids]
    as.numeric(-expm1(-exp(linpred)))
}
#' @method residuals nbhd_fit
#' @rdname nbhd_fit
#' @export
residuals.nbhd_fit = function(object, ...) {
    as.numeric(object$y) - fitted.nbhd_fit(object)
}

#' Plot coefficient estimates
#'
#' 50% and 90% credible intervals plotted by default.
#'
#' @param x a `nbhd_fit` object from [neighborhood_model()].
#' @param y ignored
#' @param intercept if `FALSE`, don't plot the intercept estimate.
#' @param inner_prob the inner credible interval probability
#' @param outer_prob the inner credible interval probability
#' @param ... Ignored.
#'
#' @return A ggplot.
#'
#' @method plot nbhd_fit
#' @export
plot.nbhd_fit = function(x, y=NULL, intercept=FALSE,
                         inner_prob = 0.5, outer_prob = 0.9, ...) {
    rlang::check_installed("ggplot2")
    probs = c((1-inner_prob)/2, 1-(1-inner_prob)/2,
              (1-outer_prob)/2, 1-(1-outer_prob)/2)
    sum_d = summary(x$post[c("coefs", "alpha")],
                    q = ~ quantile2(., probs=probs),
                    med = median) %>%
        mutate(variable = if_else(str_starts(variable, "coefs"),
                                  str_sub(variable, 7, -2),
                                  variable))
    if (isFALSE(intercept))
        sum_d = filter(sum_d, variable != "(Intercept)")
    names(sum_d)[2:5] = c("in_l", "in_h", "out_l", "out_h")

    ggplot2::ggplot(sum_d, ggplot2::aes(variable, med)) +
        ggplot2::geom_hline(yintercept=0, lty="dashed", color="#444444") +
        ggplot2::geom_linerange(ggplot2::aes(ymin=out_l, ymax=out_h), size=0.75) +
        ggplot2::geom_linerange(ggplot2::aes(ymin=in_l, ymax=in_h), size=1.25) +
        ggplot2::geom_point(size=2) +
        ggplot2::coord_flip() +
        ggplot2::labs(x=NULL, y="Estimate")
}

#' @rdname nbhd_fit
#' @export
as.matrix.nbhd_fit = function(x, ...) {
    posterior::as_draws_matrix(nm$post)
}
#' @rdname nbhd_fit
#' @export
as.data.frame.nbhd_fit = function(x, ...) {
    posterior::as_draws_df(nm$post)
}
