.register = function() {
    rstantools::rstan_config(".")
    cpp11::cpp_register(".")
    cpp11 = readr::read_lines("src/cpp11.cpp", lazy=FALSE)
    insert_header = c("#ifdef RCPP_USE_GLOBAL_ROSTREAM",
                      "Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();",
                      "Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();",
                      "#endif", "",
                      "RcppExport SEXP _rcpp_module_boot_stan_fit4nbhd_glmm_mod();", "")
    insert_call =  '    {"_rcpp_module_boot_stan_fit4nbhd_glmm_mod",  (DL_FUNC) &_rcpp_module_boot_stan_fit4nbhd_glmm_mod,  0},'
    header_idx = which(stringr::str_detect(cpp11, '#include "cpp11/declarations.hpp"'))[1] + 1
    call_idx = which(stringr::str_detect(cpp11, '\\s+\\{NULL, NULL, 0\\}'))[1] - 1
    force_idx = which(stringr::str_detect(cpp11, "R_forceSymbols"))[1]
    new_cpp11 = c(cpp11[1:header_idx],
                  insert_header,
                  cpp11[(header_idx+1):call_idx],
                  insert_call,
                  cpp11[(call_idx+1):(force_idx-1)],
                  cpp11[(force_idx+1):length(cpp11)])
    readr::write_lines(new_cpp11, "src/cpp11.cpp")
}

.compile = function(clean=F) {
    .register()
    if (clean) pkgbuild::clean_dll()
    pkgbuild::compile_dll(compile_attributes=FALSE)
    devtools::load_all(".")
}
