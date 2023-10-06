#' @export
M_PER_MI = 1609.34


#' Filter block data to a radius of a FIPS code
#'
#' @param fips the FIPS code to center the area at
#' @param block_d the block data, with a `$fips` column matching `fips` argument.
#' @param geom_d the block geometry data
#' @param distance the radius of the area, in miles
local_area = function(fips, block_d, geom_d, dist=0.5) {
    dist = dist * M_PER_MI
    idx = match(fips, block_d$fips)
    stopifnot(!is.na(idx))
    dists = as.numeric(s2_distance_matrix(block_d$centroid[idx],
                                          block_d$centroid))
    filter(block_d, dists <= dist) %>%
        left_join(geom_d, by=intersect(names(.), names(geom_d))) %>%
        st_as_sf(sf_column_name="geometry")
}



#' Binned Residual Plot
#'
#' @param model the model object, which should have `fitted`, `resid`, etc. methods.
#'
#' @export
binned_resid = function(model, bins=16) {
    y_fit = fitted(model)
    y_res = resid(model, type="response")
    tibble(fitted = cut(y_fit, bins),
           y_res = y_res) %>%
        group_by(fitted) %>%
        summarize(mean = mean(y_res),
                  n = n()) %>%
        mutate(fitted = parse_number(as.character(fitted))) %>%
    ggplot(aes(x=fitted, y=mean, size=n)) +
        geom_point()
}


# Helpful function adapted from: https://stat.ethz.ch/pipermail/r-help/2005-September/079872.html
#' Calculate AUC
#'
#' @param x the predictor
#' @param y a binary indicator
#'
#' @returns the scalar AUC
fastAUC <- function(x, y) {
    x1 = x[y==1]; n1 = length(x1);
    x2 = x[y==0]; n2 = length(x2);
    r = rank(c(x1, x2))
    (sum(r[1:n1]) - n1*(n1 + 1)/2) / (n1 * n2)
}
