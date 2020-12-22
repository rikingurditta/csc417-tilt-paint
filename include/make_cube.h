#include <Eigen/Core>

//Despite name, actually makes triangle mesh of any axis-aligned rectangular prism
//Input:
//  corner - position of corner with smallest x, y, z coordinates
//  size - length, height, width of rectangular prism
//Output:
//  V - 8 by 3 positions of vertices of rectangular prism
//  F - 12 by 3 faces of rectangular prism
void make_cube(const Eigen::Vector3d &corner, const Eigen::Vector3d &size, Eigen::MatrixXd &V, Eigen::MatrixXi &F);