#include "naive_viscosity_2d.h"

void naive_viscosity_2d(const Eigen::VectorXd &u, const int nx, const int ny, Eigen::VectorXd &u_new) {
    int vel_dx_grid_size = (nx + 1) * ny;
    int vel_dy_grid_size = nx * (ny + 1);
    u_new = Eigen::VectorXd::Zero(vel_dx_grid_size + vel_dy_grid_size);
    for (int i = 1; i < nx; i++) {
        for (int j = 1; j < ny - 1; j++) {
            u_new(i + j * (nx + 1)) = (u(i - 1 + (j - 1) * (nx + 1))
                                       + u(i + (j - 1) * (nx + 1))
                                       + u(i + 1 + (j - 1) * (nx + 1))
                                       + u(i - 1 + j * (nx + 1))
                                       + u(i + j * (nx + 1))
                                       + u(i + 1 + j * (nx + 1))
                                       + u(i - 1 + (j + 1) * (nx + 1))
                                       + u(i + (j + 1) * (nx + 1))
                                       + u(i + 1 + (j + 1) * (nx + 1))) / 9;
        }
    }
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny; j++) {
            u_new(vel_dx_grid_size + i + j * nx) = (u(vel_dx_grid_size + i - 1 + (j - 1) * nx)
                                                    + u(vel_dx_grid_size + i + (j - 1) * nx)
                                                    + u(vel_dx_grid_size + i + 1 + (j - 1) * nx)
                                                    + u(vel_dx_grid_size + i - 1 + j * nx)
                                                    + u(vel_dx_grid_size + i + j * nx)
                                                    + u(vel_dx_grid_size + i + 1 + j * nx)
                                                    + u(vel_dx_grid_size + i - 1 + (j + 1) * nx)
                                                    + u(vel_dx_grid_size + i + (j + 1) * nx)
                                                    + u(vel_dx_grid_size + i + 1 + (j + 1) * nx)) / 9;
        }
    }
}