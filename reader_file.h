#ifndef READER_FILE_H
#define READER_FILE_H

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include "quadtree_t.h"
#include <dlib/image_io.h>

static void load_input(hough_settings &settings, quadtree_t &node) {

    node.m_level = 0;
    node.m_root = &node;

    dlib::load_png(node.im, settings.input_filename);

    node.sat = SAT(node.im.nr(), node.im.nc());

    node.m_indexes.reserve(node.im.nr() * node.im.nc());


    for (int j = 0, index = 0; j < node.im.nr(); j++) {
        for (int i = 0; i < node.im.nc(); i++) {

            node.sat.SAT_X(j, i) = node.sat.getValueAt(node.sat.SAT_X, j, i - 1) + node.sat.getValueAt(node.sat.SAT_X, j - 1, i) - node.sat.getValueAt(node.sat.SAT_X, j - 1, i - 1);
            node.sat.SAT_Y(j, i) = node.sat.getValueAt(node.sat.SAT_Y, j, i - 1) + node.sat.getValueAt(node.sat.SAT_Y, j - 1, i) - node.sat.getValueAt(node.sat.SAT_Y, j - 1, i - 1);
            node.sat.SAT_Z(j, i) = node.sat.getValueAt(node.sat.SAT_Z, j, i - 1) + node.sat.getValueAt(node.sat.SAT_Z, j - 1, i) - node.sat.getValueAt(node.sat.SAT_Z, j - 1, i - 1);
            node.sat.SAT_XX(j, i) = node.sat.getValueAt(node.sat.SAT_XX, j, i - 1) + node.sat.getValueAt(node.sat.SAT_XX, j - 1, i) - node.sat.getValueAt(node.sat.SAT_XX, j - 1, i - 1);
            node.sat.SAT_XY(j, i) = node.sat.getValueAt(node.sat.SAT_XY, j, i - 1) + node.sat.getValueAt(node.sat.SAT_XY, j - 1, i) - node.sat.getValueAt(node.sat.SAT_XY, j - 1, i - 1);
            node.sat.SAT_XZ(j, i) = node.sat.getValueAt(node.sat.SAT_XZ, j, i - 1) + node.sat.getValueAt(node.sat.SAT_XZ, j - 1, i) - node.sat.getValueAt(node.sat.SAT_XZ, j - 1, i - 1);
            node.sat.SAT_YY(j, i) = node.sat.getValueAt(node.sat.SAT_YY, j, i - 1) + node.sat.getValueAt(node.sat.SAT_YY, j - 1, i) - node.sat.getValueAt(node.sat.SAT_YY, j - 1, i - 1);
            node.sat.SAT_YZ(j, i) = node.sat.getValueAt(node.sat.SAT_YZ, j, i - 1) + node.sat.getValueAt(node.sat.SAT_YZ, j - 1, i) - node.sat.getValueAt(node.sat.SAT_YZ, j - 1, i - 1);
            node.sat.SAT_ZZ(j, i) = node.sat.getValueAt(node.sat.SAT_ZZ, j, i - 1) + node.sat.getValueAt(node.sat.SAT_ZZ, j - 1, i) - node.sat.getValueAt(node.sat.SAT_ZZ, j - 1, i - 1);

            node.sat.SAT_SAMPLES(j, i) = node.sat.getValueAt(node.sat.SAT_SAMPLES, j, i - 1) + node.sat.getValueAt(node.sat.SAT_SAMPLES, j - 1, i) - node.sat.getValueAt(node.sat.SAT_SAMPLES, j - 1, i - 1);

            double z = node.im[j][i];
            if (z > 0) {
                double x = (i - settings.camera_cx) * z * settings.inv_camera_fx;
                double y = (j - settings.camera_cy) * z * settings.inv_camera_fy;

                node.sat.SAT_X(j, i) += x;
                node.sat.SAT_Y(j, i) += y;
                node.sat.SAT_Z(j, i) += z;

                node.sat.SAT_XX(j, i) += x * x;
                node.sat.SAT_YY(j, i) += y * y;
                node.sat.SAT_ZZ(j, i) += z * z;

                node.sat.SAT_XY(j, i) += x * y;
                node.sat.SAT_XZ(j, i) += x * z;
                node.sat.SAT_YZ(j, i) += y * z;

                node.sat.SAT_SAMPLES(j, i) += 1;

                node.m_indexes.push_back(index);
                node.m_pts.push_back(dvec3(x, y, z));
                settings.max_point_distance = std::max(settings.max_point_distance, dlib::length(node.m_pts[index]));
                node.m_ij.push_back(ivec2(i,j));
                node.m_colors.push_back(dvec3(0, 0, 0));
                index++;
            }
        }
    }

    node.m_middle(0) = node.sat.SAT_SAMPLES.nc() >> 1;
    node.m_middle(1) = node.sat.SAT_SAMPLES.nr() >> 1;

    node.top_left_bounds(0) = 0;
    node.top_left_bounds(1) = 0;
    node.bot_right_bounds(0) = node.sat.SAT_SAMPLES.nc() - 1;
    node.bot_right_bounds(1) = node.sat.SAT_SAMPLES.nr() - 1;

}

#endif // READER_FILE_H
