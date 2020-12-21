#include "advect_2d.h"

void advect_2d(double dt, Eigen::MatrixXd &particles, Eigen::MatrixXd &particle_velocities) {
    particles += particle_velocities * dt;
}