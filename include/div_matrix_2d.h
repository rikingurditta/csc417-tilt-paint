#include <Eigen/Core>
#include <Eigen/Sparse>

//Input:
//  nx - grid x size
//  ny - grid y size
//Output:
//  B - sparse matrix so that if g = (∂f/∂x ∂f/∂y)^T for f defined on a grid, then 1/Δx * Bg = ∇·f
void div_matrix_2d(int nx, int ny, Eigen::SparseMatrix<double> B);