#include "include/paint_headers.h"
#include <iostream>

void make_grid(
    const int nx,
    const int ny,
    const Eigen::RowVector3d corner,
    const double spacing,
    Eigen::MatrixXd & grid_points)
{
    // Shoutout to Alec Jacobson's CSC419 course for giving us staggered grid code:
    //https://github.com/alecjacobson/geometry-processing-mesh-reconstruction/blob/master/src/poisson_surface_reconstruction.cpp

    // padding: number of cells beyond bounding box of input points
    const double pad = 8;

    // Compute positions of grid nodes
    grid_points.resize(nx*ny, 3);

        for(int i = 0; i < nx; i++) {
            for(int j = 0; j < ny; j++){
                    // Convert subscript to index
                    int ind = nx * i + j;
                    grid_points.row(ind) = corner + spacing * Eigen::RowVector3d(i,j, 0.0);
            }
        }

}