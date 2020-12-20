#include "include/paint_headers.h"
#include <iostream>
#include <Eigen/Core>
#include <cstdlib>


#include "pressure_projection_2d.h"
#include "fixed_zero_velocities_2d.h"
#include "div_matrix_2d.h"
#include "grad_matrix_2d.h"


int main(int argc, char *argv[]) {
    // Main architecture goes as follows:
    // Values are hardcoded below, unless otherwise specified by command line arguments
    Eigen::RowVector3d corner = Eigen::RowVector3d(0.1, 0.1, 0.1);
    int nx = 20, ny = 50;
    double spacing = 1.0;



    // Define a grid by the origin, dimensions in x, y,
    // a corner, and the grid spacing
    Eigen::MatrixXd grid_points;
    make_grid(nx, ny, corner, spacing, grid_points);

    // print block
    // for (int i = 0; i < 5; i++) {
    //     std::cout << grid_points.row(i)<< std::endl;
    // }

    std::cout << "Hello, World!" << std::endl;

    int vel_dx_grid_size = (nx + 1) * ny;
    int vel_dy_grid_size = nx * (ny + 1);
    Eigen::VectorXd u = Eigen::VectorXd::Zero(vel_dx_grid_size + vel_dy_grid_size);
    Eigen::VectorXd u_new = Eigen::VectorXd::Zero(vel_dx_grid_size + vel_dy_grid_size);
    // create matrix B which, given gradient, calculates divergence
    Eigen::SparseMatrix<double> B;
    div_matrix_2d(nx, ny, B);
    // create matrix PP which gets rid of boundary velocities
    Eigen::SparseMatrix<double> PP;
    fixed_zero_velocities(nx, ny, PP);
    // create gradient matrix D
    Eigen::SparseMatrix<double> D;
    grad_matrix_2d(nx, ny, D);
    pressure_projection_2d(u, nx, ny, spacing, 1., PP, B, D, u_new);
    return EXIT_SUCCESS;
}
