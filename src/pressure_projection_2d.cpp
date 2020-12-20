#include <iostream>
#include "pressure_projection_2d.h"

void pressure_projection_2d(const Eigen::VectorXd &u,
                            const int nx,
                            const int ny,
                            const double dx,
                            const double rho,
                            const Eigen::SparseMatrix<double> &PP,
                            const Eigen::SparseMatrix<double> &B,
                            const Eigen::SparseMatrix<double> &D,
                            Eigen::VectorXd &u_new) {
    // solve linear system D * B * p = ρ/dt² B * PP * u, i.e. ∇·∇p = ∇·u
    Eigen::SparseMatrix<double> A = B * D;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    std::cout << "compute success: " << (solver.info() == Eigen::Success) << "\n";
    Eigen::VectorXd p = solver.solve(B * PP * u);
    p *= rho / (dx * dx);
    std::cout << "solve success: " << (solver.info() == Eigen::Success) << "\n";
}