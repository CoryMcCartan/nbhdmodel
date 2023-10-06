#include <queue>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <cpp11.hpp>

using namespace cpp11;
namespace writable = cpp11::writable;

/*
 * Subsets to the relevant region for the map model
 *
 * Returns a list with 3 elements:
 *      idx, the block indices for the relevant blocks
 *      incl, the neighborhood inclusion indicator
 *      rings, a vector of 'rings', i.e., graph distances from the starting block
 */
[[cpp11::register]]
list find_relevant(integers nbhd, list graph) {
    std::unordered_set<int> in_nbhd(nbhd.begin(), nbhd.end());
    std::unordered_set<int> checked;
    checked.insert(nbhd[0]);


    std::vector<int> idx;
    idx.reserve(nbhd.size());
    std::vector<bool> incl;
    incl.reserve(nbhd.size());
    std::vector<int> rings;
    rings.reserve(nbhd.size());

    std::queue<int> to_check;
    to_check.push(nbhd[0]);

    // breadth-first search
    int level = 0;
    while (!to_check.empty()) {
        int n_level = to_check.size();

        while (n_level--) { // for all nodes in this level/ring
            int vtx = to_check.front();
            to_check.pop();

            integers nbors = graph[vtx - 1];
            bool any = false;
            for (int nbor : nbors) {
                // if this node has a neighbor in the neighborhood
                if (!any && in_nbhd.find(nbor) != in_nbhd.end()){
                    any = true;
                    break;
                }
            }

            if (any) {
                idx.push_back(vtx);
                rings.push_back(level);
                incl.push_back(in_nbhd.find(vtx) != in_nbhd.end());

                // check neighbors
                for (int nbor : nbors) {
                    // if not already checked, add to queue
                    if (checked.find(nbor) == checked.end()) {
                        to_check.push(nbor);
                        checked.insert(nbor);
                    }
                }

            }
        }

        level++;
    }

    return writable::list({"idx"_nm = idx,
                          "incl"_nm = incl,
                          "rings"_nm = rings});
}

/*
 * Having calculated distances, check again which nodes should be included
 * under the model, and filter out those that should not be
 *
 * Returns a doubles vector that shows the fraction of nearer neighbors that are
 * in the neighborhood. This will be positive if the node should be kept
 */
[[cpp11::register]]
doubles check_incl(integers idx, logicals incl, integers rings,
                    doubles dists, list graph) {
    int n = idx.size();
    writable::doubles out(n);

    std::unordered_map<int, int> lookup;
    for (int i = 0; i < n; i++) {
        lookup[idx[i]] = i;
    }

    for (int i = 0; i < n; i++) {
        integers nbors = graph[idx[i] - 1];
        int ring = rings[i];
        double dist = dists[i];

        double ok = 0;
        int nearer = 0;
        for (int nbor : nbors) {
            auto j_it =  lookup.find(nbor);
            if (j_it == lookup.end()) continue; // nbor is outside check area
            int j = j_it->second;

            if (rings[j] < ring || (rings[j] == ring && dists[j] <= dist)) {
                nearer++;
                ok += (double) (int) incl[j];
            }
        }

        out[i] = ok / nearer;
    }

    return out;
}
