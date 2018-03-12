#include "hough.h"
#include "peak_detection.h"
#include "voting.h"
#include "plane_t.h"
#include <vector>
#include <iostream>
#include <chrono>
#include "accumulatorball_t.h"
#include <logger.h>
#include <dlib/image_io.h>

accumulatorball_t *kht3d( std::vector<plane_t> &planes, quadtree_t &father, hough_settings &settings, vector<kernel_t> &used_kernels) {

    Logger log_full = Logger("Hough Transform").start_counting();

    Logger log_step;

   /**
     *
     * Subdividing Procedure
     *
     * */

    log_step = Logger("Subdivision").start_counting();
    father.subdivide(settings);
    log_step.end_counting();

    /**
      *
      * Initializes the Accumulator
      *
      * */

    accumulatorball_t *accum = new accumulatorball_t(settings.max_point_distance, settings.n_rho, settings.n_phi);

    /**
      *
      * Voting Procedure
      *
      * */

    log_step = Logger("Voting").start_counting();
    std::vector<bin_t> used_bins;
    voting(father, *accum, used_bins, used_kernels);
    log_step.end_counting();

    /**
      *
      * Peak Detection Procedure
      *
      * */

    log_step = Logger("Peak Detection").start_counting();
    peak_detection(planes, *accum, used_kernels, used_bins);
    log_step.end_counting();

    log_full.end_counting();

   /**
     *
     * Sorting planes by representativeness
     *
     * */

    std::sort(planes.begin(), planes.end());

    /**
     *
     * Coloring planes and points
     *
     * */

    for (unsigned int i = 0; i < planes.size(); i++) {
        dvec3 cor;
        switch(i % 6) {
            case 0:
                cor = dvec3((int)(255/(int)(i/6+1)), 0, 0)/255.0;
                break;
            case 1:
                cor = dvec3(0, (int)(255/(int)(i/6+1)), 0)/255.0;
                break;
            case 2:
                cor = dvec3(0, 0, (int)(255/(int)(i/6+1)))/255.0;
                break;
            case 3:
                cor = dvec3(0, (int)(255/(int)(i/6+1)), (int)(255/(int)(i/6+1)))/255.0;
                break;
            case 4:
                cor = dvec3((int)(255/(int)(i/6+1)), 0, (int)(255/(int)(i/6+1)))/255.0;
                break;
            case 5:
                cor = dvec3((int)(255/(int)(i/6+1)), (int)(255/(int)(i/6+1)), 0)/255.0;
                break;
        }
        planes[i].m_color = cor;

        for (quadtree_t *node : planes[i].nodes) {
            node->color = cor;
        }
    }

    int half_window = 1;

    std::vector<dvec3> normals(father.m_indexes.size());
    for (int idx : father.m_indexes) {
        if (father.m_pts[idx](2) > 0) {
            int x_min = father.m_ij[idx](0) - half_window;
            int y_min = father.m_ij[idx](1) - half_window;
            int x_max = father.m_ij[idx](0) + half_window+1;
            int y_max = father.m_ij[idx](1) + half_window+1;

            quadtree_t *v_node = new quadtree_t();
            v_node->m_root = &father;

            v_node->top_left_bounds = ivec2(x_min, y_min);
            v_node->bot_right_bounds = ivec2(x_max, y_max);

            v_node->compute_samples();

            if (v_node->samples >= 3 ) {
                v_node->compute_centroid();
                v_node->run_PCA();
                if (v_node->normal(2) > 0) {
                    v_node->normal *= -1;
                }
                normals[idx] = v_node->normal;
            }
        }
    }

    dlib::array2d<dlib::rgb_pixel> output(father.sat.SAT_SAMPLES.nr(), father.sat.SAT_SAMPLES.nc());
    for (int j = 0; j < output.nr(); j++) {
        for (int i = 0; i < output.nc(); i++) {
            output[j][i] = dlib::rgb_pixel(0, 0, 0);
        }
    }

    for (plane_t plane : planes) {
        dvec3 cor = plane.m_color * 255;
        for (int idx : father.m_indexes) {

            if (plane.distance2plane(father.m_pts[idx]) < settings.max_distance2plane ) {
                if (plane.m_normal(2) > 0) {
                    plane.m_normal *= -1;
                }
                if (dlib::dot(plane.m_normal, normals[idx]) > 0.65) {
                    output[father.m_ij[idx](1)][father.m_ij[idx](0)] = dlib::rgb_pixel(cor(0), cor(1), cor(2));
                }
            }
        }
    }

    dlib::save_png(output, settings.output_filename);
    return accum;
}
