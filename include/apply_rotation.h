#include <Eigen/Core>

void apply_rotation(const Eigen::MatrixXd &R, const Eigen::RowVector3d &centre, const Eigen::MatrixXd &points,
                    Eigen::MatrixXd &result);