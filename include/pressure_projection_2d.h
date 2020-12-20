#include <Eigen/Core>
#include <Eigen/Sparse>

//Input:
//  u - current velocity grid as a stacked vector (x_vel y_vel)^T where x_vel is the horizontal staggered velocity grid
//      of size nx+1 by ny and y_vel is is the horizontal staggered velocity grid of size nx by ny+1
//  nx - grid x size
//  ny - grid y size
//  PP - sparse matrix so that PP * u = u without edge velocities
//  B - sparse divergence matrix
//  D - sparse gradient matrix
//Output:
//  u_new - new velocity grid
void pressure_projection_2d(const Eigen::VectorXd &u,
                            const int nx,
                            const int ny,
                            const double dx,
                            const double rho,
                            const Eigen::SparseMatrix<double> &PP,
                            const Eigen::SparseMatrix<double> &B,
                            const Eigen::SparseMatrix<double> &D,
                            Eigen::VectorXd &u_new);