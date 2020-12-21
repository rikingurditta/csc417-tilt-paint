#include "apply_gravity_2d.h"

void apply_gravity_2d(double dt, const Eigen::Vector2d g, Eigen::MatrixXd &particle_velocities) {
    particle_velocities += dt * g.transpose().replicate(particle_velocities.rows(), 1);
}