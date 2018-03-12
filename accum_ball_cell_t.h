#ifndef ACCUM_BALL_CELL_T
#define ACCUM_BALL_CELL_T

#include "accum_cell_t.h"

// An accumulator cell (θ, φ) with an array of (ρ) cells
class accum_ball_cell_t {
    public:

    accum_ball_cell_t(int size) {
        bins = new accum_cell_t[size]();
    }

    ~accum_ball_cell_t() {
        delete [] bins;
    }

    accum_cell_t *bins;

private:
    accum_ball_cell_t(const accum_ball_cell_t &);
};

#endif
