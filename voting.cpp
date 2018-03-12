#include "voting.h"
#include <ctime>
#include <limits>

inline bool cast_vote(std::vector<bin_t> &used_bins, accum_cell_t &cell, kernel_t &kernel, double votes, double t, int p, int r) {

    // Cluster representativeness
    votes = votes * kernel.node->representativeness;

    // Test if the cell has already been voted by this kernel
    if (cell.verify_cell(kernel.node)) {
        // Test if the previous vote was smaller than the current
        if (cell.last_casted_vote < votes) {
            // Remove
            cell.bin += -cell.last_casted_vote+votes;
            cell.last_casted_vote = votes;
        } else {
            return false;
        }
    // First node vote
    } else {
        // Store how many votes will be cast
        cell.last_casted_vote = votes;
        // Increment votes
        cell.bin += votes;
        // Add reference of the node that votes for this cell
        cell.add_reference(kernel.node);

        if (!cell.voted) {
            used_bins.push_back(bin_t(t, p, r));
            cell.voted = true;
        }
        // Track last node that votes for this cell
        cell.apply_cell(kernel.node);
    }
    return true;
}

// Voting in 1 dimensions (rho)
inline bool gaussian_vote_1d(accumulatorball_t &accum, kernel_t &kernel, std::vector<bin_t> &used_bins, 
                             const double theta_index, const int phi_index, const int rho_start_index,
                             const double theta, const double phi) {
   
    int rho_index, p, inc_rho_index;
    double rho, t, inc_rho, gauss;
    bool voted = false;

    double votes;

    p = phi_index;
    t = theta_index;
    inc_rho = accum.m_delta_rho;
    inc_rho_index = 1;


    // Voting in the RHO array the direction "positive"
    for (rho_index = rho_start_index, rho = 0.0 ;; rho_index += inc_rho_index, rho += inc_rho) {
        if (rho_index < 0) {
            inc_rho_index *= -1;
        }
        if (accum.process_rho(t, p, rho_index) == false) {
            break;
        }
        gauss = kernel.trivariated_gaussian_dist_normal(rho, phi, theta);

        if (gauss < kernel.voting_limit) {
            break;
        }
        if ((votes = gauss) >= 0.0) {
            voted = cast_vote(used_bins, accum.at(t, p, rho_index), kernel, votes, t, p, rho_index);
        } else {
            break;
        }
    }


    // Voting in the RHO array the direction "negative"
    inc_rho_index = -1;
    for (rho_index = rho_start_index-1, rho = -inc_rho;;rho_index += inc_rho_index, rho -= inc_rho) {
        if (rho_index < 0) {
            inc_rho_index *= -1;
        }
        if (accum.process_rho(t, p, rho_index) == false) {
            break;
        }

        gauss = kernel.trivariated_gaussian_dist_normal(rho, phi, theta);
        if (gauss < kernel.voting_limit) {
            break;
        }

        if ((votes = gauss) >= 0.0) {
            voted = cast_vote(used_bins, accum.at(t, p, rho_index), kernel, votes, t, p, rho_index);
        } else {
            break;
        }
    }
    return voted;
}

// Voting in 2 dimensions (theta, phi)
void gaussian_vote_2d( accumulatorball_t &accum, kernel_t &kernel, std::vector<bin_t> &used_bins, 
                       const double theta_start_index, int phi_start_index, int rho_start_index,
                       const double phi_start, int inc_phi_index) {
    double inc_phi = (double)inc_phi_index * accum.m_delta_angle;
    int phi_index, rho_index = rho_start_index;
    bool voted, voting_ended = false;
    double theta_index, theta, phi, theta_init;

    // Voting in the PHI array in the direction "inc_phi_index"
    for (phi_index = phi_start_index, phi = phi_start; !voting_ended ; phi_index += inc_phi_index, phi += inc_phi) {

        theta_index = theta_start_index;
        if (accum.process_phi(theta_index, phi_index)) {
            inc_phi_index *= -1;
        }
        theta_init = accum.fix_theta(theta_index,phi_index);

        double inc_theta = accum.delta_theta(phi_start_index);
        double inc_theta_index = accum.delta_theta_index(phi_index);
        double curr_theta_index;
        voted = false;
        voting_ended = true;

        // Voting in the THETA array in the direction "positive"
        for (curr_theta_index = theta_init, theta = 0.0 ;; curr_theta_index+=inc_theta_index, theta += inc_theta) {
            accum.process_theta(curr_theta_index);
            accum.initialize(curr_theta_index, phi_index);
            voted = gaussian_vote_1d(accum, kernel, used_bins, curr_theta_index, phi_index, rho_index, theta, phi);
            if (voted) {
                voting_ended = false;
            } else {
                break;
            }
        }

       // Voting in the THETA array in the direction "negative"
        for (curr_theta_index = theta_init - inc_theta_index, theta = -inc_theta;; curr_theta_index -= inc_theta_index, theta -= inc_theta) {
            accum.process_theta(curr_theta_index);
            accum.initialize(curr_theta_index, phi_index);
            voted = gaussian_vote_1d(accum, kernel, used_bins, curr_theta_index, phi_index, rho_index, theta, phi);
            if (voted) {
                voting_ended = false;
            } else{
                break;
            }
       }
    }
}

// Voting in 3 dimensions (theta, phi and rho)
void gaussian_vote_3d(kernel_t &kernel, accumulatorball_t &accum, std::vector<bin_t> &used_bins) {
    accum.at(kernel.theta_index, kernel.phi_index, kernel.rho_index).top = true;

    gaussian_vote_2d(accum, kernel, used_bins, kernel.theta_index, kernel.phi_index, kernel.rho_index, 0, +1);

    int phi_index = kernel.phi_index-1;
    double theta_index = kernel.theta_index;
    accum.process_phi(theta_index, phi_index);

    gaussian_vote_2d(accum, kernel, used_bins, theta_index, phi_index, kernel.rho_index, -accum.m_delta_angle, -1);
    }

// Calculates gaussian kernel parameters
void kernel_calculation(quadtree_t &node, accumulatorball_t &accum, std::vector<kernel_t> &used_kernels) {
    kernel_t kernel;
    kernel.node = &node;

    if (dlib::dot(dlib::normalize(node.m_centroid), node.normal) < 0) {
        node.normal *= -1.0;
    }
    kernel.phi = acos(node.normal(2));
    kernel.theta = atan2(node.normal(1), node.normal(0));
    kernel.rho = dlib::dot(node.m_centroid, node.normal);
    accum.process_limits(kernel.theta_index, kernel.phi_index, kernel.rho_index);
    accum.get_index(kernel.theta, kernel.phi, kernel.rho,  kernel.theta_index, kernel.phi_index, kernel.rho_index);
    kernel.theta_index = accum.fix_theta(kernel.theta_index, kernel.phi_index);

    kernel.kernel_load_parameters();
    used_kernels.push_back(kernel);
}

void voting(quadtree_t &root, accumulatorball_t &accumulator, std::vector<bin_t> &used_bins, std::vector<kernel_t> &used_kernels) {
    int i = 0;
    for (quadtree_t *node : root.coplanarNodes) {
        kernel_calculation(*node, accumulator, used_kernels);
        gaussian_vote_3d(used_kernels[i++], accumulator, used_bins);
    }
}




