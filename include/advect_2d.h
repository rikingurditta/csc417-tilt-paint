#include <Eigen/Core>

//Advects particles, i.e. moves them
//Input:
//  dt - time step
//  nx - grid x size
//  ny - grid y size
//  dx - grid cell width
//  particle-velocities - n by 2 velocities of each particle
//  particles - n by 2 current positions of each particle
//Output:
//  particles - n by 2 updated positions of each particle
void advect_2d(double dt, int nx, int ny, double dx, Eigen::MatrixXd &particle_velocities, Eigen::MatrixXd &particles);