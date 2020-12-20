#include "include/paint_headers.h"
#include <iostream>

 // Helper function for turning grid nodes to indices
int grid_to_i(
  const int i, 
  const int j, 
  const int k, 
  const int n_x, 
  const int n_y) {
    return i + n_x * (j + k * n_y);
};

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////

  // std::cout <<" nx "<< nx <<std::endl;    
  // std::cout <<" corner "<< corner <<std::endl;    



  std::vector<Eigen::Triplet<double>> tripletList;
  // tripletList.reserve(nx*ny*nz);
  W.resize(P.rows(), nx*ny*nz);

  // Now for each point, calculate the corresponding weight for each component
  // If corner is the origin, get the point w.r.t grid

  // Just iterate over each row
  for (int i = 0; i < P.rows(); i++) {
    Eigen::RowVector3d h_vector, p_row, point_wrt_grid, snap_point, percent_denom, snap_percent, 
    floor_grid_comp, weight_row;
    h_vector = Eigen::RowVector3d(h, h, h);
    p_row = P.row(i);
    point_wrt_grid = p_row - corner;   

    // Find "lowest" grid point to snap to. Divide by h
    snap_point = point_wrt_grid.cwiseProduct(h_vector.cwiseInverse());
    // std::cout << "B"<< snap_point <<std::endl;    

    floor_grid_comp =  Eigen::RowVector3d(floor(snap_point[0]), floor(snap_point[1]), floor(snap_point[2]));
    snap_percent = snap_point - floor_grid_comp;

    // Get percentages and snap points
    int spx = floor_grid_comp[0];
    int spy = floor_grid_comp[1];
    int spz = floor_grid_comp[2];
  
    // std::cout << spx <<" "<< spy<<" "<< spz <<std::endl;    

    double px = snap_percent[0];
    double py = snap_percent[1];
    double pz = snap_percent[2];

    // Construct an array to get the index of each grid node, and another array to keep track of the 
    // percentages for that grid node.
    int box_i[] = {
      grid_to_i(spx, spy, spz, nx, ny), // 0

      grid_to_i(spx + 1, spy, spz, nx, ny), // 1
      grid_to_i(spx, spy + 1, spz, nx, ny), // 2
      grid_to_i(spx, spy, spz + 1, nx, ny), // 3

      grid_to_i(spx + 1, spy + 1, spz, nx, ny), // 4
      grid_to_i(spx + 1, spy, spz + 1, nx, ny), // 5
      grid_to_i(spx, spy + 1, spz + 1, nx, ny), // 6

      grid_to_i(spx + 1, spy + 1, spz + 1, nx, ny), //7
    };

    double percent_i[] = {
      (1.0 - px)*(1.0 - py)*(1.0 - pz), // 0

      px*(1.0 - py)*(1.0 - pz), // 1
      (1.0 - px)*py*(1.0 - pz), // 2
      (1.0 - px)*(1.0 - py)*pz, // 3

      px*py*(1.0 - pz), // 4
      px*(1.0 - py)*pz, // 5
      (1.0 - px)*py*pz, // 6

      px*py*pz, //7
    };

    // Iterate over all 8 points and push them into the triplet list
    for (int n = 0; n < 8; n ++) {
      double first = i;
      double second = box_i[n];
      double third = percent_i[n];
      // std::cout << first <<" "<< second<<" "<< third <<std::endl;    
      tripletList.push_back(Eigen::Triplet<double>(first, second, third));
      
      // Try coeffref
      // std::cout << first <<" "<< second<<" "<< third <<std::endl;    
      // W.coeffRef(first, second) = third;
    }
  }
  W.setFromTriplets(tripletList.begin(), tripletList.end());

  // std::cout << "W: "<< W.rows() <<"x"<< W.cols() <<std::endl;
}
