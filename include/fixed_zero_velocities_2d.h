#include <Eigen/Core>
#include <Eigen/Sparse>

//Input:
//  nx - grid x size
//  ny - grid y size
//Output:
//  PP - sparse matrix that zeroes out velocity at edges of the grid
//       so if f = (u v)^T where u is horizontal velocity grid and v is vertical velocity grid, then PP * f returns
//       u and v stacked, but with the edges replaced with 0s
void fixed_zero_velocities(int nx, int ny, Eigen::SparseMatrix<double> &PP);