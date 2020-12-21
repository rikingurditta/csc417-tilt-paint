#include <Eigen/Core>

// despite name, can make any axis-aligned rectangular prism
void make_cube(const Eigen::Vector3d &corner, const Eigen::Vector3d &size, Eigen::MatrixXd &V, Eigen::MatrixXi &F);