#include "include/paint_headers.h"
#include <iostream>
#include <igl/cat.h>


void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
  

  Eigen::SparseMatrix<double> Dx, Dy, Dz;
  // Fill each matrix
  // Dx
  fd_partial_derivative(nx, ny, nz, h, 0, Dx);
  // std::cout << "Dx size: "<< Dx.rows()<<"x"<<Dx.cols() <<std::endl;
  // Dy
  fd_partial_derivative(nx, ny, nz, h, 1, Dy);
  // std::cout << "Dy size: "<< Dy.rows()<<"x"<<Dy.cols() <<std::endl;
  // Dz
  fd_partial_derivative(nx, ny, nz, h, 2, Dz);
  // std::cout << "Dz size: "<< Dz.rows()<<"x"<<Dz.cols() <<std::endl;

  //concatenate Dx, Dy, Dz to fill G
  Eigen::SparseMatrix<double> DxDy;
  // Put Dx and DY into DxDY
  igl::cat(1, Dx, Dy, DxDy);
  // Put temporary DxDy and Dz together, put into G
  igl::cat(1, DxDy, Dz, G);

  // std::cout << "G size: "<< G.rows()<<"x"<<G.cols() <<std::endl;

}
