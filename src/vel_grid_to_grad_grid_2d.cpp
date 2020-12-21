#include "vel_grid_to_grad_grid_2d.h"
#include <iostream>

void vel_grid_to_grad_grid_2d(int nx, int ny, Eigen::SparseMatrix<double> &PP) {
    std::vector<Eigen::Triplet<double>> tl;
    int in_dx_grid_size = (nx + 1) * ny;
    int in_dy_grid_size = nx * (ny + 1);
    int out_dx_grid_size = (nx - 1) * ny;
    int out_dy_grid_size = nx * (ny - 1);
    std::cout << out_dx_grid_size + out_dy_grid_size << ", " << in_dx_grid_size + in_dy_grid_size << "\n";
    tl.reserve(out_dx_grid_size + out_dy_grid_size);
    PP.reserve(out_dx_grid_size + out_dy_grid_size);
    PP.resize(out_dx_grid_size + out_dy_grid_size, in_dx_grid_size + in_dy_grid_size);
    for (int j = 0; j < ny; j++) {
        for (int i = 1; i < nx; i++) {  // skips over i=0 and i=nx+1 because we don't want edges of grid
            tl.emplace_back((i - 1) + j * (nx - 1), i + j * (nx + 1), 1);
            // PP.coeffRef((i - 1) + j * (nx - 1), i + j * (nx + 1)) = 1;
            // std::cout << (i - 1) + j * (nx - 1) << ", " << i + j * (nx + 1) << "\n";
        }
    }
    for (int j = 1; j < ny; j++) {  // skips over j=0 and i=ny+1 because we don't want edges of grid
        for (int i = 0; i < nx; i++) {
            tl.emplace_back(out_dx_grid_size + i + (j - 1) * nx, in_dx_grid_size + i + j * nx, 1);
            // PP.coeffRef(out_dx_grid_size + i + (j - 1) * nx, in_dx_grid_size + i + j * nx) = 1;
            // std::cout << out_dx_grid_size + i + (j - 1) * nx << ", " << in_dx_grid_size + i + j * nx << "\n";
        }
    }
    PP.setFromTriplets(tl.begin(), tl.end());
}