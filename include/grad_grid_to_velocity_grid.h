#include <Eigen/Core>
#include <Eigen/Sparse>

//Input:
//  nx - grid x size
//  ny - grid y size
//Output:
//  D_to_vel - sparse matrix so that if df is a gradient grid, then D_to_vel * df is the same grid but padded so that it
//             has the same size as a velocity grid
void grad_grid_to_velocity_grid(int nx, int ny, Eigen::SparseMatrix<double> &B);