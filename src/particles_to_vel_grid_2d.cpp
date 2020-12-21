#include "particles_to_vel_grid_2d.h"
#include <iostream>

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
        int xgrid_x = floor(x);
        int xgrid_y = floor(y - 0.5);
        int ygrid_x = floor(x - 0.5);
        int ygrid_y = floor(y);
//        std::cout << grid_x << ", " << grid_y << "\n";
        // bilinear interpolation weights
        double x_w00 = (xgrid_x + 1 - x) * (xgrid_y + 1 - y);
        double x_w10 = (x - xgrid_x) * (xgrid_y + 1 - y);
        double x_w01 = (xgrid_x + 1 - x) * (y - xgrid_y);
        double x_w11 = (x - xgrid_x) * (y - xgrid_y);
        if (0 <= xgrid_y and xgrid_y < ny) {
            if (0 <= xgrid_x and xgrid_x < nx + 1) {
                u(xgrid_x + xgrid_y * (nx + 1)) += particle_velocities(i, 0) * x_w00 * dx * dx / particle_volume;
                u(xgrid_x + (xgrid_y + 1) * (nx + 1)) += particle_velocities(i, 0) * x_w01 * dx * dx / particle_volume;
            }
            if (-1 <= xgrid_x and xgrid_x < nx) {
                u(xgrid_x + 1 + xgrid_y * (nx + 1)) += particle_velocities(i, 0) * x_w10 * dx * dx / particle_volume;
                u(xgrid_x + 1 + (xgrid_y + 1) * (nx + 1)) += particle_velocities(i, 0) * x_w11 * dx * dx / particle_volume;
            }
        }
        double y_w00 = (ygrid_x + 1 - x) * (ygrid_y + 1 - y);
        double y_w10 = (x - ygrid_x) * (ygrid_y + 1 - y);
        double y_w01 = (ygrid_x + 1 - x) * (y - ygrid_y);
        double y_w11 = (x - ygrid_x) * (y - ygrid_y);
        if (0 <= ygrid_x and ygrid_x < nx) {
            if (0 <= ygrid_y and ygrid_y < ny + 1) {
                u(vel_x_size + ygrid_x + ygrid_y * nx) += particle_velocities(i, 1) * y_w00 * dx * dx / particle_volume;
                u(vel_x_size + ygrid_x + 1 + ygrid_y * nx) += particle_velocities(i, 1) * y_w10 * dx * dx / particle_volume;
            }
            if (-1 <= ygrid_y and ygrid_y < ny) {
                u(vel_x_size + ygrid_x + (ygrid_y + 1) * nx) += particle_velocities(i, 1) * y_w01 * dx * dx / particle_volume;
                u(vel_x_size + ygrid_x + 1 + (ygrid_y + 1) * nx) += particle_velocities(i, 1) * y_w11 * dx * dx / particle_volume;
            }
        }
    }
}