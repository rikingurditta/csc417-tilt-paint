#include "grad_grid_to_vel_grid_2d.h"

void grad_grid_to_vel_grid_2d(int nx, int ny, Eigen::SparseMatrix<double> &D_to_vel) {
    // grad grid is between grid cells
    int grad_x_size = (nx - 1) * ny;
    int grad_y_size = nx * (ny - 1);
    // velocity grid is around grid cells
    int vel_x_size = (nx + 1) * ny;
    int vel_y_size = nx * (ny + 1);
    D_to_vel.resize(vel_x_size + vel_y_size, grad_x_size + grad_y_size);
    std::vector <Eigen::Triplet<double>> tl;
    tl.reserve(grad_x_size + grad_y_size);
    for (int i = 0; i < nx - 1; i++) {
        for (int j = 0; j < ny; j++) {
            tl.emplace_back((i + 1) + j * (nx + 1), i + j * (nx - 1), 1);
        }
    }
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny - 1; j++) {
            tl.emplace_back(vel_x_size + i + (j + 1) * nx, grad_x_size + i + j * nx, 1);
        }
    }
    D_to_vel.setFromTriplets(tl.begin(), tl.end());
}