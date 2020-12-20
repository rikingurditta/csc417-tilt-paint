#include "make_cube.h"

void make_cube(const Eigen::Vector3d &corner, const Eigen::Vector3d &size, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
    V = Eigen::MatrixXd::Zero(8, 3);
    F = Eigen::MatrixXi::Zero(12, 3);
    Eigen::Vector3d sx = Eigen::Vector3d::Zero();
    sx(0) = size(0);
    Eigen::Vector3d sy = Eigen::Vector3d::Zero();
    sy(1) = size(1);
    Eigen::Vector3d sz = Eigen::Vector3d::Zero();
    sz(2) = size(2);
    V.row(0) = corner;
    V.row(1) = corner + sx;
    V.row(2) = corner + sy;
    V.row(3) = corner + sy + sx;
    V.row(4) = corner + sz;
    V.row(5) = corner + sz + sx;
    V.row(6) = corner + sz + sy;
    V.row(7) = corner + sz + sy + sx;
    // TODO: fix face orientations
    F << 0, 1, 4,
            0, 2, 6,
            0, 6, 4,
            0, 2, 3,
            0, 3, 1,
            1, 4, 5,
            1, 7, 3,
            1, 5, 7,
            4, 7, 5,
            4, 6, 7,
            3, 6, 2,
            3, 7, 6;
}