#include <Eigen/Core>
#include <Eigen/Sparse>

void interpolate_vel_2d(Eigen::MatrixXd &particles,
                        Eigen::VectorXd &u,
                        int nx, int ny,
                        double dx,
                        Eigen::MatrixXd &particle_velocities);