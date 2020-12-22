#include "advect_2d.h"

void advect_2d(double dt, int nx, int ny, double dx, Eigen::MatrixXd &particle_velocities, Eigen::MatrixXd &particles) {
    particles += particle_velocities * dt;
    // clamp particles that try to move out of grid
    for (int i = 0; i < particles.rows(); i++) {
        particles(i, 0) = fmin(nx * dx - dx / 4, fmax(dx / 4, particles(i, 0)));
        particles(i, 1) = fmin(ny * dx - dx / 4, fmax(dx / 4, particles(i, 1)));
    }
}