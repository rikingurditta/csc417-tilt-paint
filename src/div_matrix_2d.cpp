#include <div_matrix_2d.h>

void div_matrix_2d(const int nx, const int ny, const double dx, const double dy, Eigen::SparseMatrix<double> B) {
    std::vector<Eigen::Triplet<double>> tl;
    tl.reserve(4 * (nx - 1) * (ny - 1));
    // sizes of staggered grids
    int dx_grid_size = (nx - 1) * ny;
    for (int i = 0; i < nx - 1; i++) {
        for (int j = 0; j < ny - 1; j++) {
            int cell = i + j * (nx - 1);
            tl.emplace_back(cell, i + j * (nx - 1), -1 / dx);
            tl.emplace_back(cell, i + j * (nx - 1) + 1, 1 / dx);
            tl.emplace_back(cell, i + j * nx + dx_grid_size, -1 / dy);
            tl.emplace_back(cell, i + j * nx + 1 + dx_grid_size, 1 / dy);
        }
    }
    B.setFromTriplets(tl.begin(), tl.end());
}