#include "apply_gravity_2d.h"

void apply_gravity_2d(double dt, const Eigen::Vector2d g, Eigen::MatrixXd &particle_velocities) {
    for (int i = 0; i < particle_velocities.rows(); i++) {
        particle_velocities.row(i) += dt * g;
    }
}