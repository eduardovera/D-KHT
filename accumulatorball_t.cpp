#include "accumulatorball_t.h"

accumulatorball_t::accumulatorball_t(const double max_distance, const int rho_num, const int phi_num) {
    neighbors_size = 27;
    m_theta_max = PI2;
    m_phi_max = PI;

    m_rho_length = rho_num;
    m_phi_length = phi_num;
    m_phi_length_half = .5 * phi_num;
    m_data.resize(m_phi_length + 1);
    m_delta_angle = PI / (double)m_phi_length;
    m_delta_rho = max_distance / (double)rho_num;
    for (int p = 0; p <= m_phi_length; p++) {
        double m_phi = (double)p / (double)m_phi_length * m_phi_max;
        int length = std::max(1, static_cast<int>(round((double)m_phi_length * 2.0 * sin(m_phi))));
        m_data[p].resize(length, NULL);
    }
}

accumulatorball_t::~accumulatorball_t() {
    for (size_t i = 0; i < m_data.size(); i++) {
        for (size_t j = 0; j < m_data[i].size(); j++) {
            if (m_data[i][j] != NULL) {
                delete m_data[i][j];
            }
        }
    }
}
