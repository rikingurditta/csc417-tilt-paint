#include <Eigen/Core>
#include <Eigen/Sparse>

//Input:
//  nx - grid x size
//  ny - grid y size
//  dx - width of a grid cell
//  dy - height of a grid cell
//Output:
//  B - sparse matrix so that if g = (∂f/∂x ∂f/∂y)^T for f defined on a grid, then Bg = ∇·x
void div_matrix_2d(int nx, int ny, double dx, double dy, Eigen::SparseMatrix<double> B);