#ifndef ACCUM_CELL_T
#define ACCUM_CELL_T


#include "quadtree_t.h"


// A cell (θ, φ, ρ) of the accumulator which stores the votes
class accum_cell_t {
    public:
        accum_cell_t() {
            last_casted_vote = 0.0;
            last_node_voted = NULL;
            visited = false;
            top = false;
            voted = false;
            peak = false;
            bin = 0.0;
        }

        ~accum_cell_t() {

        }

        inline bool verify_cell(quadtree_t * ref) {
            return (last_node_voted == ref);
        }

        inline void apply_cell(quadtree_t * ref) {
            last_node_voted = ref;
            this->bin++;
        }

        inline void add_reference(quadtree_t * ref) {
            ref_node.insert(ref);
        }

        std::set<quadtree_t *> ref_node;

        quadtree_t *last_node_voted;
        bool peak, visited, voted, top;
        double last_casted_vote, bin;
        double theta_index;
        int phi_index, rho_index;
};

#endif
