#include <queue>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <cpp11.hpp>

using namespace cpp11;
namespace writable = cpp11::writable;

/*
 * Get the precincts which are `ring` steps away from `start`
 */
[[cpp11::register]]
list get_within_ring(int ring, int start, list graph) {
    std::unordered_set<int> visited;
    visited.insert(start);

    writable::integers out_vtx(1);
    writable::integers out_ring(1);
    out_vtx[0] = start;
    out_ring[0] = 0;

    std::queue<int> to_visit;
    to_visit.push(start);
    int level = 0;
    while (level < ring && !to_visit.empty()) {
        int n_level = to_visit.size();

        while (n_level--) { // for all nodes in this level/ring
            int vtx = to_visit.front();
            to_visit.pop();

            integers nbors = graph[vtx - 1];
            for (int nbor : nbors) {
                // if not already checked, add to queue
                if (visited.find(nbor) == visited.end()) {
                    to_visit.push(nbor);
                    visited.insert(nbor);
                    out_vtx.push_back(nbor);
                    out_ring.push_back(level + 1);
                }
            }
        }

        level++;
    }

    return writable::list({"idx"_nm = out_vtx,
                          "rings"_nm = out_ring});
}

/*
 * Simulate inverse inclustion probabilities and put in a matrix
 */
[[cpp11::register]]
logicals_matrix<> sim_incl(doubles_matrix<> base_linpred, doubles ranef) {
    int nrow = base_linpred.nrow();
    int ncol = base_linpred.ncol();
    writable::logicals_matrix<> out(nrow, ncol);

    GetRNGstate();
    for (int j = 0; j < ncol; j++) {
        for (int i = 0; i < nrow; i++) {
            double linpred = base_linpred(i, j) + ranef[j];
            out(i, j) = exp_rand() >= std::exp(linpred);
        }
    }
    PutRNGstate();

    return out;
}


/*
 * Having calculated distances, check again which nodes should be included
 * under the model, and filter out those that should not be
 *
 * Returns a doubles vector that shows the fraction of nearer neighbors that are
 * in the neighborhood. This will be positive if the node should be kept
 */
[[cpp11::register]]
logicals_matrix<> fix_incl_mat(integers idx, logicals_matrix<> incl,
                               integers rings, doubles dists, list graph) {
    int n = idx.size();
    int draws = incl.ncol();

    // copy over
    writable::logicals_matrix<> out(n, draws);
    for (int j = 0; j < draws; j++) {
        for (int i = 0; i < n; i++) {
            out(i, j) = incl(i, j);
        }
    }

    std::unordered_map<int, int> lookup;
    for (int i = 0; i < n; i++) {
        lookup[idx[i]] = i;
    }

    for (int i = 0; i < n; i++) {
        integers nbors = graph[idx[i] - 1];
        int ring = rings[i];
        if (ring == 0) continue;
        double dist = dists[i];

        std::vector<bool> ok(draws, false);
        for (int nbor : nbors) {
            auto j_it =  lookup.find(nbor);
            if (j_it == lookup.end()) continue; // nbor is outside check area
            int j = j_it->second;

            if (rings[j] < ring || (rings[j] == ring && dists[j] <= dist)) {
                for (int k = 0; k < draws; k++) {
                    if (out(j, k) == TRUE) ok[k] = true;
                }
            }
        }

        for (int j = 0; j < draws; j++) {
            if (!ok[j]) {
                out(i, j) = false;
            }
        }
    }

    return out;
}
