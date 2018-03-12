#ifndef VOTING_H
#define VOTING_H

#include "accumulatorball_t.h"
#include "kernel_t.h"
#include "bin_t.h"


void voting(quadtree_t &root, accumulatorball_t &accum, std::vector<bin_t> &used_bins, std::vector<kernel_t> &used_kernels);

#endif
