#include "particles_to_vel_grid_2d.h"

void particles_to_vel_grid_2d(Eigen::MatrixXd &particles,
                              Eigen::MatrixXd &particle_velocities,
                              int nx, int ny,
                              double dx, double particle_volume,
                              Eigen::VectorXd &u) {
    int vel_x_size = (nx + 1) * ny;
    int vel_y_size = nx * (ny + 1);
    u = Eigen::VectorXd::Zero(vel_x_size + vel_y_size);
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
        // distribute velocity to grid with bilinear weights
        if (0 <= grid_x <= nx and 0 <= grid_y < ny - 1) {
            u(grid_x + grid_y * (nx + 1)) += particle_velocities(i, 0) * w00 * dx * dx / particle_volume;
            u(grid_x + 1 + grid_y * (nx + 1)) += particle_velocities(i, 0) * w10 * dx * dx / particle_volume;
            u(grid_x + (grid_y + 1) * (nx + 1)) += particle_velocities(i, 0) * w01 * dx * dx / particle_volume;
            u(grid_x + 1 + (grid_y + 1) * (nx + 1)) += particle_velocities(i, 0) * w11 * dx * dx / particle_volume;
        }
        if (0 <= grid_x < nx - 1 and 0 <= grid_y <= ny) {
            u(grid_x + grid_y * nx) += particle_velocities(i, 1) * w00 * dx * dx / particle_volume;
            u(grid_x + 1 + grid_y * nx) += particle_velocities(i, 1) * w10 * dx * dx / particle_volume;
            u(grid_x + (grid_y + 1) * nx) += particle_velocities(i, 1) * w01 * dx * dx / particle_volume;
            u(grid_x + 1 + (grid_y + 1) * nx) += particle_velocities(i, 1) * w11 * dx * dx / particle_volume;
        }
    }
}