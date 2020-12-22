#include <Eigen/Core>

//Updates particle velocities using given gravity vector
//Input:
//  dt - time step
//  g - 2d acceleration due to gravity with respect to the simulation grid
//  particle-velocities - n by 2 current velocities of each particle
//Output:
//  particle-velocities - n by 2 updated velocities of each particle
void apply_gravity_2d(double dt, Eigen::Vector2d g, Eigen::MatrixXd &particle_velocities);