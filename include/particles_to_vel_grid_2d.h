#include <Eigen/Core>

void particles_to_vel_grid_2d(Eigen::MatrixXd &particles,
                                  Eigen::MatrixXd &particle_velocities,
                                  int nx, int ny,
                                  double dx, double particle_volume,
                                  Eigen::VectorXd &u);