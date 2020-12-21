#include <Eigen/Core>

//naive, not physics based, just turns each velocity into average of velocities around it
//this is problematic as it causes the viscosity to depend very heavily on timestep
//Input:
//  u - current velocity grid as a stacked vector (x_vel y_vel)^T where x_vel is the horizontal staggered velocity grid
//      of size nx+1 by ny and y_vel is is the horizontal staggered velocity grid of size nx by ny+1
//  nx - grid x size
//  ny - grid y size
//Output:
//  u_new - new velocity grid with viscosity applied
void naive_viscosity_2d(const Eigen::VectorXd &u, const int nx, const int ny, Eigen::VectorXd &u_new);