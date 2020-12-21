#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/OrderingMethods>

//Input:
//  u - current velocity grid as a stacked vector (x_vel y_vel)^T where x_vel is the horizontal staggered velocity grid
//      of size nx+1 by ny and y_vel is is the horizontal staggered velocity grid of size nx by ny+1
//  nx - grid x size
//  ny - grid y size
//  nx - grid cell width
//  dt - time step
//  rho - fluid density
//  PP - sparse matrix so that PP * u = u without edge velocities
//  B - sparse divergence matrix
//  D - sparse gradient matrix
//Output:
//  u_new - new velocity grid
void pressure_projection_2d(const Eigen::VectorXd &u,
                            const int nx,
                            const int ny,
                            const double dx,
                            const double dt,
                            const double rho,
                            const Eigen::SparseMatrix<double> &PP,
                            const Eigen::SparseMatrix<double> &B,
                            const Eigen::SparseMatrix<double> &D,
                            const Eigen::SparseMatrix<double> &D_to_vel,
                            Eigen::VectorXd &u_new);