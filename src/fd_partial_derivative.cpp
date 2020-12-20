#include "include/paint_headers.h"
#include <iostream>

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
  
  // std::cout << "nx: "<< nx <<std::endl;
  // std::cout << "ny: "<< ny <<std::endl;
  // std::cout << "nz: "<< nz <<std::endl;
  // std::cout << "h: "<< h <<std::endl;

  // Check to see dir's value, if dir = 0 output Dx, 1 output Dy, 2 output Dz

  // Dx
  if (dir == 0) {
    D.conservativeResize((nx-1)*ny*nz, (nx*ny*nz));
    // std::cout << "D size: "<< D.rows()<<"x"<<D.cols() <<std::endl;

    // Now fill in grid
    for(int i = 0; i < nx-1; i++){
      for(int j = 0; j < ny; j++){
        for(int k = 0; k < nz; k++){
          // You are now at a staggered grid coordinate here
          const auto ind_one = i + (nx - 1)*(j + k*ny);
          const auto ind_two = i + nx * (j + k * ny);
          const auto ind_three = (i + 1) + nx * (j + k*ny);
          D.coeffRef(ind_one, ind_two) = -1.0;
          D.coeffRef(ind_one, ind_three) = 1.0;
        }
      }
    }
  }

  // Dy
  if (dir == 1) {
    D.conservativeResize((nx*(ny-1)*nz), nx*ny*nz);
    // Now fill in grid
    for(int i = 0; i < nx; i++){
      for(int j = 0; j < ny-1; j++){
        for(int k = 0; k < nz; k++){
          // You are now at a staggered grid coordinate here
          const auto ind_one = i + nx*(j + k*(ny - 1));
          const auto ind_two = i + nx * (j + k * ny);
          const auto ind_three = i + nx * ((j + 1) + k*ny);
          D.coeffRef(ind_one, ind_two) = -1.0;
          D.coeffRef(ind_one, ind_three) = 1.0;
        }
      }
    }
  }

// Dz
  if (dir == 2) {
    D.conservativeResize((nx*ny*(nz-1)), nx*ny*nz);
    // Now fill in grid
    for(int i = 0; i < nx; i++){
      for(int j = 0; j < ny; j++){
        for(int k = 0; k < nz-1; k++){
          // You are now at a staggered grid coordinate here
          const auto ind_one = i + nx * (j + k * ny);
          const auto ind_two = i + nx * (j + k * ny);
          const auto ind_three = i + nx * (j + (k + 1) *ny);
          D.coeffRef(ind_one, ind_two) = -1.0;
          D.coeffRef(ind_one, ind_three) = 1.0;
        }
      }
    }
  }
}
