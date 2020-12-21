#include "grad_matrix_2d.h"

void grad_matrix_2d(const int nx, const int ny, Eigen::SparseMatrix<double> &D) {
    std::vector<Eigen::Triplet<double>> tl;
    tl.reserve(2 * (nx - 1) * ny + nx * (ny - 1));
    int dx_grid_size = (nx - 1) * ny;
    int dy_grid_size = nx * (ny - 1);
    for (int x = 0; x < nx - 1; x++) {
        for (int y = 0; y < ny; y++) {
            // x on staggered grid is between x and x+1 on main grid
            tl.emplace_back(x + y * (nx - 1), x + y * nx, -1);
            tl.emplace_back(x + y * (nx - 1), x + 1 + y * nx, 1);
        }
    }
    for (int x = 0; x < nx; x++) {
        for (int y = 0; y < ny - 1; y++) {
            // y on staggered grid is between y and y+1 on main grid
            tl.emplace_back(dx_grid_size + x + y * nx, x + y * nx, -1);
            tl.emplace_back(dx_grid_size + x + y * nx, x + (y + 1) * nx, 1);
        }
    }
    D.resize(dx_grid_size + dy_grid_size, nx * ny);
    D.setFromTriplets(tl.begin(), tl.end());
}