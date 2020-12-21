#include <Eigen/Core>

void advect_2d(double dt, int nx, int ny, double dx, Eigen::MatrixXd &particles, Eigen::MatrixXd &particle_velocities);