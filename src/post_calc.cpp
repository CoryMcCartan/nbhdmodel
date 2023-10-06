#include <vector>
#include <cmath>
#include <cpp11.hpp>

using namespace cpp11;
namespace writable = cpp11::writable;
/*
 * Simulate inverse inclustion probabilities and put in a matrix
 */
[[cpp11::register]]
doubles pr_incl(doubles_matrix<> base_linpred, doubles ranef) {
    int nrow = base_linpred.nrow();
    int ncol = base_linpred.ncol();
    writable::doubles out(nrow);

    for (int i = 0; i < nrow; i++) {
        double accuml = 0;
        for (int j = 0; j < ncol; j++) {
            double linpred = base_linpred(i, j) + ranef[j];
            accuml += std::exp(-std::exp(linpred));
        }
        out[i] = accuml / (double) ncol;
    }

    return out;
}
