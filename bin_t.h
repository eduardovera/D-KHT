#ifndef BIN_T
#define BIN_T

// A structure that stores the coordinates of an accumulator cell
struct bin_t {
    bin_t() { }
    bin_t(double theta, short phi, short rho) : theta_index(theta), phi_index(phi), rho_index(rho) { }

    double theta_index;
    short phi_index;
    short rho_index;
    double votes;

    bool operator < (const bin_t bin) const {
        return (votes > bin.votes);
    }
};

#endif
