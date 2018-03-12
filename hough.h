#ifndef HOUGH_H
#define HOUGH_H

#include "accumulatorball_t.h"
#include "kernel_t.h"
#include "plane_t.h"

#include <vector>

class hough_settings;
class quadtree_t;

accumulatorball_t *kht3d(std::vector<plane_t> &planes, quadtree_t &father, hough_settings &settings, vector<kernel_t> &used_kernels);

#endif
