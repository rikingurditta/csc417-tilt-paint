#include "fixed_zero_velocities.h"

void fixed_zero_velocities(int nx, int ny, Eigen::SparseMatrix<double> PTP) {
    std::vector<Eigen::Triplet<double>> tl;
    int dx_grid_size = (nx + 1) * ny;
    int dy_grid_size = nx * (ny + 1);
    tl.reserve(dx_grid_size + dy_grid_size - 2 * (nx + ny));
    for (int i = 1; i < nx; i++) {  // skips over i=0 and i=nx+1 because sides of grid should be 0
        for (int j = 0; j < ny; j++) {
            tl.emplace_back(i + j * (nx + 1), i + j * (nx + 1), 1);
        }
    }
    for (int i = 0; i < nx; i++) {
        for (int j = 1; j < ny; j++) {  // skips over j=0 and i=ny+1 because top and bottom of grid should be 0
            tl.emplace_back(dx_grid_size + i + j * nx, dx_grid_size + i + j * nx, 1);
        }
    }
    PTP.resize(dx_grid_size + dy_grid_size, dx_grid_size + dy_grid_size);
    PTP.setFromTriplets(tl.begin(), tl.end());
}