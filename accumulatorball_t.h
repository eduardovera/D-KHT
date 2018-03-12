#ifndef ACCUMULATORBALL_T
#define ACCUMULATORBALL_T

#include "accum_ball_cell_t.h"
#include "bin_t.h"
#include <cstring>
#include <vector>


class accumulatorball_t {
    public:
        accumulatorball_t(const double max_distance, const int rho_num, const int phi_num);
        ~accumulatorball_t();

        inline bool visited_neighbour(std::set<accum_cell_t*> &neighbors) {
            for (accum_cell_t * cell : neighbors) {
                if (cell->visited) {
                    return true;
                }
            }
            return false;
        }

        inline std::set<quadtree_t *> convolution_nodes(const double theta_index, const short phi_index, const short rho_index) {
            std::set<quadtree_t *> nodes;
            std::set<accum_cell_t*> neighbors = get_neighbors(theta_index, phi_index, rho_index, 26);
            for (accum_cell_t * cell : neighbors) {
                std::move(cell->ref_node.begin(), cell->ref_node.end(), std::inserter(nodes, nodes.end()));
            }
            return nodes;
        }

        inline void set_visited(const double theta_index, const short phi_index, const short rho_index) {
            std::set<accum_cell_t*> neighbors = get_neighbors(theta_index, phi_index, rho_index, 27);
            for (accum_cell_t * cell : neighbors) {
                cell->visited = true;
            }
        }
   
        inline void initialize(double theta_index, int phi_index) {
            int t = get_theta_index(theta_index,phi_index);
            if (m_data[phi_index][t] == NULL) {
                m_data[phi_index][t] =  new accum_ball_cell_t(m_rho_length);
            }
        }

        inline double convolution_value(const double theta_index, const int phi_index, const int rho_index) {
            double acc_value = 0.0;
            std::set<accum_cell_t*> neighbors = get_neighbors(theta_index, phi_index, rho_index, 6);
            accum_cell_t *center = &at(theta_index, phi_index, rho_index);
            acc_value = center->bin * 0.2002;
            for (accum_cell_t *n : neighbors ) {
                acc_value += n->bin * 0.1333;
            }
            return acc_value;
        }

        inline void set_visited(std::vector<accum_cell_t*> &neighbors) {
            for (accum_cell_t * cell : neighbors) {
                cell->visited = true;
            }
        }

        inline double delta_theta(const double &phi_index) {
            return m_theta_max / (double)(m_data[phi_index].size());
        }

        inline double delta_theta_index(const double &phi_index) {
            return 1.0 / (double)(m_data[phi_index].size());
        }

        inline void process_theta(double &theta_index) {
            theta_index = (theta_index+1.0) - (int)(theta_index+1.0);
        }

        inline bool process_phi(double &theta_index, int &phi_index) {
            process_theta(theta_index);
            if (phi_index < 0) {
                phi_index = std::abs(phi_index);
                theta_index = (theta_index + .5) - (int)(theta_index + .5);
                return true;
            } else if (phi_index > m_phi_length) {
                phi_index = m_phi_length + m_phi_length - phi_index;
                theta_index = (theta_index + .5) - (int)(theta_index + .5);
                return true;
            }
            return false;
        }

        inline bool process_rho(double &theta_index, int &phi_index, int &rho_index) {
            if (rho_index < 0) {
                rho_index = std::abs(rho_index);
                phi_index = m_phi_length - phi_index;
                theta_index = (theta_index+.5) - (int)(theta_index+.5);
                return true;
            } else if (rho_index >= m_rho_length) {
                return false;
            }
            return true;
        }

        inline bool visited_neighbour(const double theta_index, const short phi_index, const short rho_index)
        {
           std::set<accum_cell_t*> neighbors = get_neighbors(theta_index,phi_index,rho_index, 26);
           for (accum_cell_t * cell : neighbors) {
              if (cell->visited) return true;
           }
           return false;
        }

        inline std::set<accum_cell_t*> get_neighbors( const double theta_index, const short phi_index, const short rho_index, const int neighborhood_size) {
            std::set<accum_cell_t*> result;

            int p, r;
            double t, theta;

            if (phi_index == 0 || phi_index == m_phi_length) {
                theta = 0.5;
            } else {
                theta = theta_index;
            }

            //                                 center[1]   /  direct-linked[7]  /         /          semi-direct-linked[19]          /      /     diagonal-linked[27]    /
            static const short offset_x[26] = {          0,  0,  0,  0, +1, -1,      +1, -1, +1, -1, +1, -1, -1, +1,  0,  0,  0,  0,     +1, +1, +1, -1, +1, -1, -1, -1,};
            static const short offset_y[26] = {         +1, -1,  0,  0,  0,  0,      +1, -1, -1, +1,  0,  0,  0,  0, +1, -1, -1, +1,     +1, +1, -1, +1, -1, -1, +1, -1,};
            static const short offset_z[26] = {          0,  0, +1, -1,  0,  0,       0,  0,  0,  0, +1, -1, +1, -1, +1, -1, +1, -1,     +1, -1, +1, +1, -1, +1, -1, -1,};

            for (short i = 1; i != neighborhood_size; ++i) {
                t = theta;
                p = phi_index + offset_y[i];
                r = rho_index + offset_z[i];
                process_phi(t, p);
                t = fix_theta(t, p);
                t += (delta_theta_index(p) * (double)((offset_x[i])));

                if (!process_limits(t, p, r)) {
                    continue;
                }

                accum_cell_t *cell = &at(t, p, r);
                result.insert(cell);
            }
            return result;
        }

        inline bool process_limits(double &theta_index, int &phi_index, int &rho_index) {
            process_phi(theta_index, phi_index);
            return process_rho(theta_index, phi_index, rho_index);
        }

        inline accum_cell_t & at (const double theta, const short phi , const short rho) {
            int t = get_theta_index(theta, phi);
            if (m_data[phi][t] == NULL) {
                m_data[phi][t] = new accum_ball_cell_t(m_rho_length);
            }
            return m_data[phi][t]->bins[rho];
        }

        inline int get_theta_index(const double &theta, const int &phi_index) {
            return static_cast<int>(round(theta * (double)(m_data[phi_index].size()))) % m_data[phi_index].size();
        }

        inline void get_index(const double &theta, const double &phi, const double &rho, double &theta_index, int &phi_index, int &rho_index) {
            theta_index = (theta/PI2) + 0.5;
            phi_index = round(phi / m_delta_angle);
            rho_index = round(rho / m_delta_rho);
        }

        inline void get_values(double &theta, double &phi, double &rho, const double &theta_index, const int &phi_index, const int &rho_index) {
            theta = (theta_index-0.5) * PI2;
            phi =   (double)(phi_index) * m_delta_angle;
            rho =   (double)(rho_index) * m_delta_rho;
        }

        inline void spherical_to_cartesian(dvec3 &normal, const double &theta, const double &phi, const double &rho) {
            normal(0) = sin(phi) * cos(theta) * rho;
            normal(1) = sin(phi) * sin(theta) * rho;
            normal(2) = cos(phi) * rho;
        }

        inline double fix_theta(double theta, const int phi) {
            double p_size = m_data[phi].size();
            int t = round((theta) * p_size);
            return (t == 1) ? (0.0) : t/p_size;
        }

        int neighbors_size;

        double m_theta_max;
        double m_phi_max;

        short m_phi_length;         // Accumulator size (phi dimension)
        short m_phi_length_half;
        short m_rho_length;         // Accumulator size (rho dimension)
        short m_theta_length;       // Accumulator size (theta dimension)

        double m_delta_rho;         // Discretization step to the distance (rho)
        double m_delta_angle;       // Discretization step to the angles (theta & phi)

        std::vector< std::vector <accum_ball_cell_t*> > m_data;
};

#endif
