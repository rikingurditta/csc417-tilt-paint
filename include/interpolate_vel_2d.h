#include <Eigen/Core>
#include <Eigen/Sparse>

//Uses bilinear interpolation to distribute velocities from grid to particles
//Input:
//  particles - n by 2 matrix of positions of particles
//  u - current velocity grid as a stacked vector (x_vel y_vel)^T where x_vel is the horizontal staggered velocity grid
//      of size nx+1 by ny and y_vel is is the horizontal staggered velocity grid of size nx by ny+1
//  nx - grid x size
//  ny - grid y size
//  dx - grid cell width
//Output:
//  particle_velocities - n by 2 matrix of velocities of particles, interpolated from u
void interpolate_vel_2d(Eigen::MatrixXd &particles, Eigen::VectorXd &u, int nx, int ny, double dx,
                        Eigen::MatrixXd &particle_velocities);