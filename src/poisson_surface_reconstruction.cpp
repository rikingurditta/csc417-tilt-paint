#include "include/paint_headers.h"
#include <igl/copyleft/marching_cubes.h>
#include <algorithm>
#include <iostream>

void poisson_surface_reconstruction(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & N,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F)
{
  ////////////////////////////////////////////////////////////////////////////
  // Construct FD grid, CONGRATULATIONS! You get this for free!
  ////////////////////////////////////////////////////////////////////////////
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
  // Compute positions of grid nodes
  Eigen::MatrixXd x(nx*ny*nz, 3);
  for(int i = 0; i < nx; i++) 
  {
    for(int j = 0; j < ny; j++)
    {
      for(int k = 0; k < nz; k++)
      {
         // Convert subscript to index
         const auto ind = i + nx*(j + k * ny);
         x.row(ind) = corner + h*Eigen::RowVector3d(i,j,k);
      }
    }
  }
  Eigen::VectorXd g = Eigen::VectorXd::Zero(nx*ny*nz);

  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////

  // TODO REMOVE FUNCTION CALLS
  // Eigen::SparseMatrix<double> W;
  // fd_interpolate(nx, ny, nz, h, corner, P, W);
  // int dir = 2;
  // Eigen::SparseMatrix<double> D;
  // Eigen::SparseMatrix<double> G;

  // fd_partial_derivative(nx, ny, nz, h, dir, D);
  // fd_grad(nx, ny, nz, h, G);

  // SUBMIT THIS CODE
  Eigen::SparseMatrix<double> G, W, Wx, Wy, Wz;
  // Get weights in both big W matrix and wrt staggered grids
  // fd_interpolate(nx, ny, nz, h, corner, P, W);
  // Get staggered corners
  Eigen::RowVector3d c_x, c_y, c_z;
  c_x = Eigen::RowVector3d(h * 0.5, 0, 0) + corner;
  c_y = Eigen::RowVector3d(0, h * 0.5, 0) + corner;
  c_z = Eigen::RowVector3d(0, 0, h * 0.5) + corner;
  // Get weights for X Y Z
  // my set triplets from interpolate suddenly isn't working :(
  fd_interpolate(nx - 1, ny, nz, h, c_x, P, Wx);
  fd_interpolate(nx, ny - 1, nz, h, c_y, P, Wy);
  fd_interpolate(nx, ny, nz - 1, h, c_z, P, Wz);
  
  //  Project onto V
  Eigen::VectorXd Vx, Vy, Vz, V_total;
  Vx = Wx.transpose()*N.col(0);
  Vy = Wy.transpose()*N.col(1);
  Vz = Wz.transpose()*N.col(2);

  // Get G and solve
  fd_grad(nx, ny, nz, h, G);
  int x_shape = (nx - 1) * ny * nz;
  int y_shape = nx * (ny - 1) * nz;
  int z_shape = nx * ny * (nz - 1);
  V_total.resize(x_shape + y_shape + z_shape);
  V_total << Vx, Vy, Vz;

  // Time to solve
  Eigen::BiCGSTAB  <Eigen::SparseMatrix<double>> solver;
  Eigen::SparseMatrix<double> G_G = G.transpose() * G;
  Eigen::VectorXd RHS = G.transpose() * V_total;

  // std::cout<< "Made Solver"<<std::endl;    
  // std::cout<< "GG "<<G_G.rows()<<"x"<<G_G.cols()<<std::endl;
  // std::cout<< "RHS "<<RHS.rows()<<"x"<<RHS.cols()<<std::endl;



  solver.compute(G_G);
  // std::cout<< "Computed"<<std::endl;    

  g = solver.solve(RHS);
  // std::cout<< "Solved"<<std::endl;    


  // Add sigma
  Eigen::RowVectorXd points_ones = Eigen::RowVectorXd::Ones(P.rows());
  Eigen::VectorXd total_ones = Eigen::VectorXd::Ones(nx*ny*nz);

  double sigma = (1.0/P.rows()) * points_ones * W * g;
  g -= sigma * total_ones;

  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
