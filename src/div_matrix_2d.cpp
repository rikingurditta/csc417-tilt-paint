#include <div_matrix_2d.h>

void div_matrix_2d(const int nx, const int ny, const double dx, const double dy, Eigen::SparseMatrix<double> B) {
    std::vector<Eigen::Triplet<double>> tl;
    tl.reserve(4 * nx * ny);
    // sizes of staggered grids
    int dx_grid_size = (nx + 1) * ny;
    int dy_grid_size = nx * (ny + 1);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int cell = i + j * nx;
            tl.emplace_back(cell, i + j * (nx + 1), -1 / dx);
            tl.emplace_back(cell, i + j * (nx + 1) + 1, 1 / dx);
            tl.emplace_back(cell, i + j * nx + dx_grid_size, -1 / dy);
            tl.emplace_back(cell, i + (j + 1) * nx + dx_grid_size, 1 / dy);
        }
    }
    B.resize(nx * ny, dx_grid_size + dy_grid_size);
    B.setFromTriplets(tl.begin(), tl.end());
}