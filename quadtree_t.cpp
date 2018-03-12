#include "quadtree_t.h"
#include <limits>

quadtree_t::quadtree_t()  {
    m_root = NULL;
    m_children = NULL;
    coplanar = false;
    m_centroid = dvec3(0, 0, 0);
    color = dvec3(0.5, 0.5, 0.5);
    variance1 = variance2 = variance3 = 0.0;
    votes = 0;
    samples = 0;
}

void quadtree_t::clear() {
    if (m_children != NULL) {
        for (short i = 0; i < 4 ; i++) {
            m_children[i].clear();
        }
        delete [] m_children;
    }
    m_children = NULL;

    nonCoplanarNodes.clear();
    coplanarNodes.clear();
    m_indexes.clear();
    m_pts.clear();
    m_ij.clear();
    m_colors.clear();

    nonCoplanarNodes.shrink_to_fit();
    coplanarNodes.shrink_to_fit();
    m_indexes.shrink_to_fit();
    m_pts.shrink_to_fit();
    m_ij.shrink_to_fit();
    m_colors.shrink_to_fit();

    nonCoplanarNodes.~vector();
    coplanarNodes.~vector();
    m_indexes.~vector();
    m_pts.~vector();
    m_ij.~vector();
    m_colors.~vector();


}

void quadtree_t::compute_centroid() {
    SAT *sat = &(m_root->sat);

    int x_min = top_left_bounds(0);
    int y_min = top_left_bounds(1);
    int x_max = bot_right_bounds(0);
    int y_max = bot_right_bounds(1);

    m_centroid(0) = (
                        sat->getValueAt(sat->SAT_X, y_max, x_max)
                      - sat->getValueAt(sat->SAT_X, y_min, x_max)
                      - sat->getValueAt(sat->SAT_X, y_max, x_min)
                      + sat->getValueAt(sat->SAT_X, y_min, x_min)
                        );

    m_centroid(1) = (
                        sat->getValueAt(sat->SAT_Y, y_max, x_max)
                      - sat->getValueAt(sat->SAT_Y, y_min, x_max)
                      - sat->getValueAt(sat->SAT_Y, y_max, x_min)
                      + sat->getValueAt(sat->SAT_Y, y_min, x_min)
                        );

    m_centroid(2) = (
                        sat->getValueAt(sat->SAT_Z, y_max, x_max)
                      - sat->getValueAt(sat->SAT_Z, y_min, x_max)
                      - sat->getValueAt(sat->SAT_Z, y_max, x_min)
                      + sat->getValueAt(sat->SAT_Z, y_min, x_min)
                        );

    m_centroid /= samples;

}

void quadtree_t::compute_samples() {
    SAT *sat = &(m_root->sat);
    int x_min = top_left_bounds(0);
    int y_min = top_left_bounds(1);
    int x_max = bot_right_bounds(0);
    int y_max = bot_right_bounds(1);

    this->samples = (
                sat->getValueAt(sat->SAT_SAMPLES, y_max, x_max)
              - sat->getValueAt(sat->SAT_SAMPLES, y_min, x_max)
              - sat->getValueAt(sat->SAT_SAMPLES, y_max, x_min)
              + sat->getValueAt(sat->SAT_SAMPLES, y_min, x_min)
             );
}

void quadtree_t::subdivide( hough_settings &settings ) {

    this->compute_samples();

    if (samples < settings.s_ms) {
        m_root->nonCoplanarNodes.push_back(this);
        return;
    }

    run_PCA();

    double thickness = 1.96 * sqrt(variance1);
    if (thickness < settings.s_t) {
        coplanar = true;
        m_root->coplanarNodes.push_back(this);
        return;
    }

    m_children = new quadtree_t[4];

    m_children[0].top_left_bounds(0) = top_left_bounds(0);
    m_children[0].top_left_bounds(1) = top_left_bounds(1);
    m_children[0].bot_right_bounds(0) = m_middle(0);
    m_children[0].bot_right_bounds(1) = m_middle(1);


    m_children[1].top_left_bounds(0) = m_middle(0);
    m_children[1].top_left_bounds(1) = top_left_bounds(1);
    m_children[1].bot_right_bounds(0) = bot_right_bounds(0);
    m_children[1].bot_right_bounds(1) = m_middle(1);


    m_children[2].top_left_bounds(0) = top_left_bounds(0);
    m_children[2].top_left_bounds(1) = m_middle(1);
    m_children[2].bot_right_bounds(0) = m_middle(0);
    m_children[2].bot_right_bounds(1) = bot_right_bounds(1);


    m_children[3].top_left_bounds(0) = m_middle(0);
    m_children[3].top_left_bounds(1) = m_middle(1);
    m_children[3].bot_right_bounds(0) = bot_right_bounds(0);
    m_children[3].bot_right_bounds(1) = bot_right_bounds(1);

    int child_size_w = (bot_right_bounds(0) - top_left_bounds(0) + 1) >> 1;
    int child_size_h = (bot_right_bounds(1) - top_left_bounds(1) + 1) >> 1;

    int half_sizeX = child_size_w >> 1;
    int half_sizeY = child_size_h >> 1;

    m_children[0].m_middle(0) = m_middle(0) - half_sizeX;
    m_children[1].m_middle(0) = m_middle(0) + half_sizeX;
    m_children[2].m_middle(0) = m_middle(0) - half_sizeX;
    m_children[3].m_middle(0) = m_middle(0) + half_sizeX;

    m_children[0].m_middle(1) = m_middle(1) - half_sizeY;
    m_children[1].m_middle(1) = m_middle(1) - half_sizeY;
    m_children[2].m_middle(1) = m_middle(1) + half_sizeY;
    m_children[3].m_middle(1) = m_middle(1) + half_sizeY;

    for (int i = 0; i < 4 ; i++) {
        m_children[i].m_level = m_level + 1;
        m_children[i].m_root = m_root;
        m_children[i].subdivide(settings);
    }
}

dlib::matrix<double, 3, 3> quadtree_t::compute_covariance_matrix() {
    SAT *sat = &(m_root->sat);

    dlib::matrix<double, 3, 3> covariance(3, 3);

    int x_min = top_left_bounds(0);
    int y_min = top_left_bounds(1);
    int x_max = bot_right_bounds(0);
    int y_max = bot_right_bounds(1);

    covariance(0, 0) = (
                       (
                          sat->getValueAt(sat->SAT_XX, y_max, x_max)
                        - sat->getValueAt(sat->SAT_XX, y_min, x_max)
                        - sat->getValueAt(sat->SAT_XX, y_max, x_min)
                        + sat->getValueAt(sat->SAT_XX, y_min, x_min)
                       ) - 2 * m_centroid(0) * (
                           sat->getValueAt(sat->SAT_X, y_max, x_max)
                         - sat->getValueAt(sat->SAT_X, y_min, x_max)
                         - sat->getValueAt(sat->SAT_X, y_max, x_min)
                         + sat->getValueAt(sat->SAT_X, y_min, x_min)
                       ) + (samples * m_centroid(0) * m_centroid(0)));


    covariance(0, 1) = (
                       (
                          sat->getValueAt(sat->SAT_XY, y_max, x_max)
                        - sat->getValueAt(sat->SAT_XY, y_min, x_max)
                        - sat->getValueAt(sat->SAT_XY, y_max, x_min)
                        + sat->getValueAt(sat->SAT_XY, y_min, x_min)
                       ) - (m_centroid(1)) * (
                           sat->getValueAt(sat->SAT_X, y_max, x_max)
                         - sat->getValueAt(sat->SAT_X, y_min, x_max)
                         - sat->getValueAt(sat->SAT_X, y_max, x_min)
                         + sat->getValueAt(sat->SAT_X, y_min, x_min)
                        )  - (m_centroid(0)) * (
                           sat->getValueAt(sat->SAT_Y, y_max, x_max)
                         - sat->getValueAt(sat->SAT_Y, y_min, x_max)
                         - sat->getValueAt(sat->SAT_Y, y_max, x_min)
                         + sat->getValueAt(sat->SAT_Y, y_min, x_min)
                        ) + (samples * m_centroid(0) * m_centroid(1)));

    covariance(0, 2) = (
                       (
                           sat->getValueAt(sat->SAT_XZ, y_max, x_max)
                         - sat->getValueAt(sat->SAT_XZ, y_min, x_max)
                         - sat->getValueAt(sat->SAT_XZ, y_max, x_min)
                         + sat->getValueAt(sat->SAT_XZ, y_min, x_min)
                        ) - (m_centroid(2)) * (
                           sat->getValueAt(sat->SAT_X, y_max, x_max)
                         - sat->getValueAt(sat->SAT_X, y_min, x_max)
                         - sat->getValueAt(sat->SAT_X, y_max, x_min)
                         + sat->getValueAt(sat->SAT_X, y_min, x_min)
                        )  - (m_centroid(0)) * (
                           sat->getValueAt(sat->SAT_Z, y_max, x_max)
                         - sat->getValueAt(sat->SAT_Z, y_min, x_max)
                         - sat->getValueAt(sat->SAT_Z, y_max, x_min)
                         + sat->getValueAt(sat->SAT_Z, y_min, x_min)
                        ) + (samples * m_centroid(0) * m_centroid(2)));

    covariance(1, 1) = (
                       (
                          sat->getValueAt(sat->SAT_YY, y_max, x_max)
                        - sat->getValueAt(sat->SAT_YY, y_min, x_max)
                        - sat->getValueAt(sat->SAT_YY, y_max, x_min)
                        + sat->getValueAt(sat->SAT_YY, y_min, x_min)
                       ) - 2 * m_centroid(1) * (
                           sat->getValueAt(sat->SAT_Y, y_max, x_max)
                         - sat->getValueAt(sat->SAT_Y, y_min, x_max)
                         - sat->getValueAt(sat->SAT_Y, y_max, x_min)
                         + sat->getValueAt(sat->SAT_Y, y_min, x_min)
                       ) + (samples * m_centroid(1) * m_centroid(1)));

    covariance(1, 2) = (
                       (
                          sat->getValueAt(sat->SAT_YZ, y_max, x_max)
                        - sat->getValueAt(sat->SAT_YZ, y_min, x_max)
                        - sat->getValueAt(sat->SAT_YZ, y_max, x_min)
                        + sat->getValueAt(sat->SAT_YZ, y_min, x_min)
                       ) - (m_centroid(1)) * (
                           sat->getValueAt(sat->SAT_Z, y_max, x_max)
                         - sat->getValueAt(sat->SAT_Z, y_min, x_max)
                         - sat->getValueAt(sat->SAT_Z, y_max, x_min)
                         + sat->getValueAt(sat->SAT_Z, y_min, x_min)
                        )  - (m_centroid(2)) * (
                           sat->getValueAt(sat->SAT_Y, y_max, x_max)
                         - sat->getValueAt(sat->SAT_Y, y_min, x_max)
                         - sat->getValueAt(sat->SAT_Y, y_max, x_min)
                         + sat->getValueAt(sat->SAT_Y, y_min, x_min)
                        ) + (samples * m_centroid(1) * m_centroid(2)));

    covariance(2, 2) = (
                       (
                          sat->getValueAt(sat->SAT_ZZ, y_max, x_max)
                        - sat->getValueAt(sat->SAT_ZZ, y_min, x_max)
                        - sat->getValueAt(sat->SAT_ZZ, y_max, x_min)
                        + sat->getValueAt(sat->SAT_ZZ, y_min, x_min)
                       ) - 2 * m_centroid(2) * (
                           sat->getValueAt(sat->SAT_Z, y_max, x_max)
                         - sat->getValueAt(sat->SAT_Z, y_min, x_max)
                         - sat->getValueAt(sat->SAT_Z, y_max, x_min)
                         + sat->getValueAt(sat->SAT_Z, y_min, x_min)
                       ) + (samples * m_centroid(2) * m_centroid(2)));


    covariance(1, 0) = covariance(0, 1);
    covariance(2, 0) = covariance(0, 2);
    covariance(2, 1) = covariance(1, 2);

    covariance /= samples;

    return covariance;
}

void quadtree_t::run_PCA() {

    compute_centroid();
    m_covariance = compute_covariance_matrix();

    dlib::eigenvalue_decomposition< dlib::matrix<double, 3, 3> > eigenvalue_decomp(m_covariance);
    dlib::matrix<double> eigenvalues_vector = eigenvalue_decomp.get_real_eigenvalues();

    int min_index = 0, max_index = 0, middle_index = 0;

    if (eigenvalues_vector(1) < eigenvalues_vector(min_index)) {
        min_index = 1;
    } else if (eigenvalues_vector(1) > eigenvalues_vector(max_index)) {
        max_index = 1;
    }

    if (eigenvalues_vector(2) < eigenvalues_vector(min_index)) {
        min_index = 2;
    } else if (eigenvalues_vector(2) > eigenvalues_vector(max_index)) {
        max_index = 2;
    }

    while (middle_index == min_index || middle_index == max_index) {
        middle_index++;
    }

    variance1 = eigenvalues_vector(min_index);
    variance2 = eigenvalues_vector(middle_index);
    variance3 = eigenvalues_vector(max_index);

    variance1 /= variance2;
    variance3 /= variance2;


    dlib::matrix<double> eigenvectors_matrix = eigenvalue_decomp.get_pseudo_v();

    normal = dvec3(eigenvectors_matrix(0, min_index),    eigenvectors_matrix(1, min_index),    eigenvectors_matrix(2, min_index));
}
