#include "include/paint_headers.h"
#include <iostream>
#include <Eigen/Core>
#include <cstdlib>


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
    return EXIT_SUCCESS;
}
