#include <Eigen/Sparse>
// Header files for painting sim


// Make regular grid 
// Inputs
//   P  #P by 3 list of input points
// Outputs:
//  grid_points : nx * ny * nz x 3 matrix of grid nodes
void make_grid(
    int nx,
    int ny,
    Eigen::Vector3d corner,
    double spacing,
    Eigen::MatrixXd & grid_points);

void paint_colours(
    const int nx,
    const int ny,
    Eigen::MatrixXd & colour_mat);