#include "advect_2d.h"

void advect_2d(double dt, Eigen::MatrixXd &particles, Eigen::MatrixXd &particle_velocities) {
    for (int i = 0; i < particles.rows(); i++) {
        particles.row(i) += particle_velocities.row(i) * dt;
    }
}