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
    Eigen::RowVector3d corner,
    double spacing,
    Eigen::MatrixXd & grid_points);

// CSC419 - A2
// Construct a gradient matrix for a finite-difference grid
//
// Inputs:
//   nx  number of grid steps along the x-direction
//   ny  number of grid steps along the y-direction
//   nz  number of grid steps along the z-direction
//   h  grid step size
// Outputs:
//   G  (nx-1)*ny*nz+ nx*(ny-1)*nz+ nx*ny*(nz-1) by nx*ny*nz sparse gradient
//     matrix: G = [Dx;Dy;Dz]
//
// See also: fd_partial_derivative.h
void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G);

// CSC419 - A2
// Construct a matrix of trilinear interpolation weights for a
// finite-difference grid at a given set of points
//
// Inputs:
//   nx  number of grid steps along the x-direction
//   ny  number of grid steps along the y-direction
//   nz  number of grid steps along the z-direction
//   h  grid step size
//   corner  list of bottom-left-front corner position of grid
//   P  n by 3 list of query point locations
// Outputs:
//   W  n by (nx*ny*nz) sparse weights matrix
//
void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W);

// CSC419 - A2
// Construct a partial derivative matrix for a finite-difference grid in a
// given direction. Derivative are computed using first-order differences onto
// a staggered grid
//
// Inputs:
//   nx  number of grid steps along the x-direction
//   ny  number of grid steps along the y-direction
//   nz  number of grid steps along the z-direction
//   h  grid step size
//   dir  index indicating direction: 0-->x, 1-->y, 2-->z
// Outputs:
//   D  m by nx*ny*nz sparse partial derivative matrix, where:
//     m = (nx-1)*ny*nz  if dir = 0
//     m = nx*(ny-1)*nz  if dir = 1
//     m = nx*ny*(nz-1)  otherwise (if dir = 2)
//
// See also: fd_partial_derivative.h
void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D);

// Takes input sample points P and input normals N 
// and gives a watertight mesh using a simplified
// version of [Kazhdan et. al 2006]
//
// Inputs:
//   P  #P by 3 list of input points
//   N  #P by 3 list of input normals associated with each point in P
// Outputs:
//   V  #V by 3 list of mesh vertex positions
//   F  #F by 3 list of mesh triangle indinces into V
//
void poisson_surface_reconstruction(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & N,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F);


//TODO: get rid of this pls
void a_function(
  const float haha);