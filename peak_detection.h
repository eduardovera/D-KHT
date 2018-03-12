#include "accumulatorball_t.h"
#include "bin_t.h"
#include "kernel_t.h"
#include "plane_t.h"
#include <cmath>
#include <vector>
#include <algorithm>

inline void peak_detection(std::vector<plane_t> &planes, accumulatorball_t &accum, std::vector<kernel_t> &used_kernels, std::vector<bin_t> &used_bins) {
        for (bin_t &bin : used_bins) {
            accum_cell_t &cell = accum.at( bin.theta_index, bin.phi_index, bin.rho_index);
            bin.votes = accum.convolution_value( bin.theta_index, bin.phi_index, bin.rho_index);

            cell.bin = bin.votes;
            cell.theta_index = bin.theta_index;
            cell.rho_index = bin.rho_index;
            cell.phi_index = bin.phi_index;
        }

        std::set<accum_cell_t*> cells;
        std::vector<quadtree_t*> nodes;
        std::set<accum_cell_t*> neighbors;

        for (kernel_t kernel : used_kernels) {
            accum_cell_t *nextBin = &accum.at( kernel.theta_index, kernel.phi_index, kernel.rho_index);
            double nextVotes = nextBin->bin;
            double currVotes = 0.0;
            do {
                accum_cell_t *currBin = nextBin;
                currVotes = currBin->bin;
                nextVotes = 0.0;

                neighbors = accum.get_neighbors(currBin->theta_index, currBin->phi_index, currBin->rho_index, 27);
                std::move(currBin->ref_node.begin(), currBin->ref_node.end(), std::inserter(nodes, nodes.end()));
                for (accum_cell_t *neighbor : neighbors) {
                    if (neighbor->bin > nextVotes) {
                        nextVotes = neighbor->bin;
                        nextBin = neighbor;
                    }
                }
            } while (nextVotes > currVotes);
            std::move(nodes.begin(), nodes.end(), std::inserter(nextBin->ref_node, nextBin->ref_node.end()));
            cells.insert(nextBin);
            nodes.clear();
            neighbors.clear();
        }

        int i = 0;
        for (accum_cell_t *cell : cells) {
            plane_t p;

            accum.get_values( p.m_theta, p.m_phi, p.m_rho, cell->theta_index, cell->phi_index, cell->rho_index);

            std::move(cell->ref_node.begin(), cell->ref_node.end(), std::inserter(p.nodes, p.nodes.end()));
            i++;
            p.representativeness = 0;

            p.votes = cell->bin;
            p.ti = cell->theta_index;
            p.pi = cell->phi_index;
            p.ri = cell->rho_index;
            p.calculate();

            planes.push_back(p);

        }
}

