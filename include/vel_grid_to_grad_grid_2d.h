#include <Eigen/Core>
#include <Eigen/Sparse>

//Opposite of grad_grid_to_vel_grid_2d
//Input:
//  nx - grid x size
//  ny - grid y size
//Output:
//  PP - sparse matrix so that if u is velocity grid with horizontal component dimensions (nx + 1) * ny and vertical
//       component dimensions nx * (ny + 1), PP * u will truncate edges, resulting in grids that have sizes
//       (nx - 1) * ny and nx * (ny - 1), just like gradient grids
void vel_grid_to_grad_grid_2d(int nx, int ny, Eigen::SparseMatrix<double> &PP);