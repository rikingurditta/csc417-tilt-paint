#include <Eigen/Core>

//Uses bilinear interpolation weights to distribute velocities from particles to grid
//Input:
//  particles - n by 2 matrix of positions of particles
//  particle_velocities - n by 2 matrix of velocities of particles, interpolated from u
//  nx - grid x size
//  ny - grid y size
//  dx - grid cell width
//  particle_volume - volume of a single particle
//Output:
//  u - current velocity grid as a stacked vector (x_vel y_vel)^T where x_vel is the horizontal staggered velocity grid
//      of size nx+1 by ny and y_vel is is the horizontal staggered velocity grid of size nx by ny+1
void particles_to_vel_grid_2d(Eigen::MatrixXd &particles,
                              Eigen::MatrixXd &particle_velocities,
                              int nx, int ny,
                              double dx, double particle_volume,
                              Eigen::VectorXd &u);