#include <div_matrix_2d.h>

void div_matrix_2d(const int nx, const int ny, Eigen::SparseMatrix<double> &B) {
    std::vector<Eigen::Triplet<double>> tl;
    tl.reserve(4 * nx * ny);
    // sizes of staggered grids
    int dx_grid_size = (nx - 1) * ny;
    int dy_grid_size = nx * (ny - 1);
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            int cell = i + j * nx;
            tl.emplace_back(cell, (i - 1) + (j - 1) * (nx - 1), -1);
            tl.emplace_back(cell, i + (j - 1) * (nx - 1), 1);
            tl.emplace_back(cell, (i - 1) + (j - 1) * nx + dx_grid_size, -1);
            tl.emplace_back(cell, (i - 1) + j * nx + dx_grid_size, 1);
        }
    }
    B.resize(nx * ny, dx_grid_size + dy_grid_size);
    B.setFromTriplets(tl.begin(), tl.end());
}