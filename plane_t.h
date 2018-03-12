#ifndef PLANE_T
#define  PLANE_T

#include <vector>
#include <settings.h>
#include <quadtree_t.h>
using namespace std;

class quadtree_t;

class plane_t {
    public:

        plane_t(void) {
            representativeness = 0.0;
            samples = 0.0;
        }
        ~plane_t(void) {
        }

        inline void calculate() {

            m_normal(0) = sin(m_phi) * cos(m_theta);
            m_normal(1) = sin(m_phi) * sin(m_theta);
            m_normal(2) = cos(m_phi);

            m_position = m_normal * m_rho;

            dvec3 c_u = dvec3(0.0, 0.0, 1.0);
            if (c_u == m_normal) {
                c_u = dvec3(1.0, 0.0, 0.0);
            }

            m_cross = dlib::normalize(dlib::pointwise_multiply(c_u, m_normal));
            m_cross2 = dlib::pointwise_multiply(m_normal, m_cross);

            m_centroid = dvec3(0.0, 0.0, 0.0);
            m_scale = dvec3(0.0 , 0.0, 0.0);
            m_showing = true;
            m_rotate = 0.0;
            for (quadtree_t *node : nodes) {
                representativeness += node->representativeness;
                m_centroid += node->m_centroid;
            }
            m_centroid /= (double) nodes.size();

        }

        inline bool operator < (const plane_t p) const {
            return (representativeness > p.representativeness);
        }

        dvec3 getClosestPointToOrigin() {
            if (m_theta >= .5 * PI || m_theta < -.5 * PI) {
                return dlib::normalize(m_normal) * -m_rho;
            }
            return dlib::normalize(m_normal) * m_rho;
        }


        double distance2plane( dvec3 &point ) {
            return std::abs(dlib::dot((point - m_position), dlib::normalize(m_normal)));
        }

        double m_theta;
        double m_phi;
        double m_rho;

        dvec3 m_cross;
        dvec3 m_cross2;
        dvec3 m_position;
        dvec3 m_centroid;
        dvec3 m_normal;
        dvec3 m_desloc;
        dvec3 m_scale;
        dvec3 m_color;

        double samples;

        std::vector<dvec3> m_points;
        std::vector<ivec2> m_ij;
        std::vector<size_t> m_indexes;


        double ti, m_rotate;
        int pi, ri;
        bool m_showing;

        std::set<quadtree_t *> nodes;
        double votes;
        double representativeness;

};

#endif
