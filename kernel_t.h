#ifndef KERNEL_T
#define KERNEL_T



#include "quadtree_t.h"
#include <math.h>

#define square(var) (var*var)


class kernel_t {
    public:

        quadtree_t *node;

        double rho;
        double theta;
        double phi;

        double constant_normal;
        double voting_limit;

        double theta_index;
        int phi_index;
        int rho_index;

        bool visited;

        int votes;

        dlib::matrix<double, 3, 3> covariance_rpt_normal;
        dlib::matrix<double, 3, 3> covariance_rpt_inv_normal;

        inline void kernel_load_parameters() {
      
            dlib::matrix<double, 3, 3> covariance_xyz = node->m_covariance;
            dlib::matrix<double, 3, 3> jacobian_normal;

            dvec3 n = normalize(node->normal);

            // Jacobian Matrix calculation
            double EPS2 = 0.00001;
            dvec3 p = n * rho;
            double w = square(p(0)) + square(p(1));
            double p2 = w + square(p(2));
            double sqrtW = sqrt(w);

            double invEPS2 = 1.0 / EPS2;

            jacobian_normal(0, 0) = n(0);
            jacobian_normal(0, 1) = n(1);
            jacobian_normal(0, 2) = n(2);
            jacobian_normal(1, 0) = (sqrtW < EPS2) ? (p(0) * p(2)) * invEPS2 : (p(0) * p(2)) / (sqrtW * p2);
            jacobian_normal(1, 1) = (sqrtW < EPS2) ?(p(1) * p(2)) * invEPS2 : (p(1) * p(2)) / (sqrtW * p2);
            jacobian_normal(1, 2) = (p2 < EPS2) ? -sqrtW * invEPS2 : (-sqrtW / p2);
            jacobian_normal(2, 0) = (w < EPS2) ? -p(1) * invEPS2 : -p(1) / w;
            jacobian_normal(2, 1) = (w < EPS2) ? p(0) * invEPS2 : p(0) / w;
            jacobian_normal(2, 2) = 0.0;

            // Uncertainty propagation
            dlib::matrix<double, 3, 3> jacobian_transposed_normal = dlib::trans(jacobian_normal);
            covariance_rpt_normal = jacobian_normal * covariance_xyz * jacobian_transposed_normal;

            // Cluster representativeness
            covariance_rpt_normal(0, 0) += NONZERO;
            covariance_rpt_inv_normal = dlib::inv(covariance_rpt_normal);
            constant_normal = root22pi32 * sqrt(std::abs(dlib::det(covariance_rpt_normal)));
            dlib::eigenvalue_decomposition< dlib::matrix<double, 3, 3> > eigenvalue_decomp(covariance_rpt_normal);
            dlib::matrix<double> eigenvalues_vector = eigenvalue_decomp.get_real_eigenvalues();
            dlib::matrix<double> eigenvectors_matrix = eigenvalue_decomp.get_pseudo_v();

            // Area importance (w_a)
            double w_a = 0.75;
            // Number-of-points importance (w_d)
            double w_d = 1 - w_a;
            node->representativeness = (double)((node->bot_right_bounds(0) - node->top_left_bounds(0)) * (node->bot_right_bounds(1) - node->top_left_bounds(1))) * w_a / (double)((node->m_root->bot_right_bounds(0) - node->m_root->top_left_bounds(0)) * (node->m_root->bot_right_bounds(1) - node->m_root->top_left_bounds(1))) +
                    (double)node->samples/(double)node->m_root->m_pts.size() * (w_d);

            int min_index = dlib::index_of_min(eigenvalues_vector);

            // Voting limit calculation (g_min)
            double n_of_standard_variations = 2.0;
            double radius = sqrt( eigenvalues_vector(min_index) ) * n_of_standard_variations;
            voting_limit = trivariated_gaussian_dist_normal( eigenvectors_matrix(0, min_index) * radius, eigenvectors_matrix(1, min_index) * radius, eigenvectors_matrix(2, min_index) * radius);

       }

        // Sampling the Trivariate Gaussian Distribution
        inline double trivariated_gaussian_dist_normal(const double rho, const double phi, const double theta) {
            dlib::vector<double, 3> displacement(rho, phi, theta);
            return (std::exp(-0.5 * (dlib::trans(displacement) * covariance_rpt_inv_normal * displacement))/constant_normal);
        }
};

#endif
