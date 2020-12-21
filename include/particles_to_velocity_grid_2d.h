#include <Eigen/Core>

void particles_to_velocity_grid_2d(Eigen::MatrixXd &particles,
                                   Eigen::MatrixXd particle_velocities,
                                   int nx, int ny,
                                   Eigen::VectorXd &u);