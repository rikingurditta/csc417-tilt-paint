

void make_grid(
    const Eigen::MatrixXd & P,
    Eigen::MatrixXd grid_points,
    Eigen::MatrixXd stagger_grid_x,
    Eigen::MatrixXd stagger_grid_y,
    Eigen::MatrixXd stagger_grid_z,
    )
{
    // Shoutout to Alec Jacobson's CSC419 course for giving us staggered grid code:
    //https://github.com/alecjacobson/geometry-processing-mesh-reconstruction/blob/master/src/poisson_surface_reconstruction.cpp
    // number of input points
    const int n = P.rows();
    // Grid dimensions
    int nx, ny, nz;
    // Maximum extent (side length of bounding box) of points
    double max_extent =
        (P.colwise().maxCoeff()-P.colwise().minCoeff()).maxCoeff();
    // padding: number of cells beyond bounding box of input points
    const double pad = 8;
    // choose grid spacing (h) so that shortest side gets 30+2*pad samples
    double h  = max_extent/double(30+2*pad);
    // Place bottom-left-front corner of grid at minimum of points minus padding
    Eigen::RowVector3d corner = P.colwise().minCoeff().array()-pad*h;
    // Grid dimensions should be at least 3 
    nx = std::max((P.col(0).maxCoeff()-P.col(0).minCoeff()+(2.*pad)*h)/h,3.);
    ny = std::max((P.col(1).maxCoeff()-P.col(1).minCoeff()+(2.*pad)*h)/h,3.);
    nz = std::max((P.col(2).maxCoeff()-P.col(2).minCoeff()+(2.*pad)*h)/h,3.);

    // Get staggered corners
    Eigen::RowVector3d c_x, c_y, c_z;
    c_x = Eigen::RowVector3d(h * 0.5, 0, 0) + corner;
    c_y = Eigen::RowVector3d(0, h * 0.5, 0) + corner;
    c_z = Eigen::RowVector3d(0, 0, h * 0.5) + corner;

    // Compute positions of grid nodes
    grid_points.resize(nx*ny*nz, 3);
    // Loop over "staggered options", x = 0, y = 1, z = 2
    int s_nx, s_ny, s_nz;
    for (int mode = 0; mode < 3; mode++) {
        if (mode == 0) {
            s_nx = nx - 1;
            s_ny = ny;
            s_nz = nz;
        }
        if (mode == 1) {
            s_nx = nx;
            s_ny = ny - 1;
            s_nz = nz;
        }
        if (mode == 2) {
            s_nx = nx;
            s_ny = ny;
            s_nz = nz - 1;
        }
        for(int i = 0; i < s_nx; i++) {
            for(int j = 0; j < s_ny; j++){
                for(int k = 0; k < s_nz; k++)
                {
                    // Convert subscript to index
                    const auto ind = i + nx*(j + k * ny);
                    grid_points.row(ind) = corner + h*Eigen::RowVector3d(i,j,k);
                    // Do staggered stuff below pls

                }
            }
        }
    }

}