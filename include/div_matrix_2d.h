#include <Eigen/Core>
#include <Eigen/Sparse>

//Input:
//  nx - grid x size
//  ny - grid y size
//Output:
//  B - sparse matrix so that if f defined on a staggered grid, then 1/Δx * Bf = ∇·f
void div_matrix_2d(int nx, int ny, Eigen::SparseMatrix<double> &B);