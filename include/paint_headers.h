// Header files for painting sim

// Make regular grid 
// Inputs
//   P  #P by 3 list of input points
// Outputs:
//  grid_points : nx * ny * nz x 3 matrix of grid nodes
void make_grid(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & N,
    Eigen::MatrixXd grid_points)
{