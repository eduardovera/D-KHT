#include <iostream>


#include "settings.h"
#include "quadtree_t.h"
#include "kernel_t.h"
#include "plane_t.h"
#include "accumulatorball_t.h"
#include "hough.h"
#include <chrono>
#include "reader_file.h"
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <limits>
#include <logger.h>

vector<plane_t> identified_planes;

hough_settings settings;
quadtree_t father;
std::vector<plane_t> planes_out;
std::vector<kernel_t> used_kernels;
accumulatorball_t *accum;


void reset() {
    planes_out.clear();
    planes_out.shrink_to_fit();

    used_kernels.clear();
    used_kernels.shrink_to_fit();

    delete accum;
    father.clear();
}

void run() {

    load_input(settings, father);

    father.compute_centroid();

    accum = kht3d(planes_out, father, settings, used_kernels);
    std::cout << planes_out.size() << " PLANES FOUND." << endl;;
}

int main(void) {
//    for (;;) {
        run();
        reset();
//    }
    return 0;
}
