#ifndef QUADTREE_T
#define QUADTREE_T

#include "dlib/matrix.h"
#include "settings.h"
#include <stack>
#include <bin_t.h>
#include <sat.h>

class quadtree_t {

    public:

        quadtree_t();

        void subdivide(hough_settings &settings);
        void run_PCA();
        void clear();
        void compute_centroid();
        void compute_samples();

        dlib::matrix<double, 3, 3> compute_covariance_matrix();

        dlib::matrix<double, 3, 3> m_covariance;
        std::vector<dvec3> m_colors, m_pts;
        std::vector<ivec2> m_ij;
        std::vector<int> m_indexes;

        dlib::array2d<uint16_t> im;

        quadtree_t * m_children, *m_root;

        dvec3 normal, m_centroid, color;
        ivec2 m_middle;
        ivec2 top_left_bounds;
        ivec2 bot_right_bounds;

        double variance1, variance2, variance3, representativeness;
        short m_level;
        bool coplanar;
        int votes;
        int samples;
        vector<quadtree_t*> coplanarNodes;
        vector<quadtree_t*> nonCoplanarNodes;

        SAT sat;

};

#endif
