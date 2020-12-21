#include "interpolate_vel_2d.h"

void interpolate_vel_2d(Eigen::MatrixXd &particles,
          Eigen::VectorXd &u,
          int nx, int ny,
          double dx,
          Eigen::MatrixXd &particle_velocities) {
    particle_velocities = Eigen::MatrixXd::Zero(particles.rows(), 2);
    int vel_x_size = (nx + 1) * ny;
    for (int i = 0; i < particles.rows(); i++) {
        double x = particles(i, 0) / dx;  // unstretch by dx
        double y = particles(i, 1) / dx;
        int grid_x = floor(x);
        int grid_y = floor(y);
        // bilinear interpolation weights
        double w00 = (grid_x + 1 - x) * (grid_y + 1 - y);
        double w10 = (x - grid_x) * (grid_y + 1 - y);
        double w01 = (grid_x + 1 - x) * (y - grid_y);
        double w11 = (x - grid_x) * (y - grid_y);
        if (0 <= grid_x and grid_x <= nx and 0 <= grid_y and grid_y < ny) {
            particle_velocities(i, 0) = w00 * u(grid_x + grid_y * (nx + 1))
                                        + w10 * u((grid_x + 1) + grid_y * (nx + 1))
                                        + w01 * u(grid_x + (grid_y + 1) * (nx + 1))
                                        + w11 * u(grid_x + (grid_y + 1) * (nx + 1));
        }
        if (0 <= grid_x and grid_x < nx and 0 <= grid_y and grid_y <= ny) {
            particle_velocities(i, 1) = w00 * u(vel_x_size + grid_x + grid_y * nx)
                                        + w10 * u(vel_x_size + (grid_x + 1) + grid_y * nx)
                                        + w01 * u(vel_x_size + grid_x + (grid_y + 1) * nx)
                                        + w11 * u(vel_x_size + grid_x + (grid_y + 1) * nx);
        }
    }
}