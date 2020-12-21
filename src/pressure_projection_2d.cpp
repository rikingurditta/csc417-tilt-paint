#include <iostream>
#include "pressure_projection_2d.h"

void pressure_projection_2d(const Eigen::VectorXd &u,
                            const int nx,
                            const int ny,
                            const double dx,
                            const double dt,
                            const double rho,
                            const Eigen::SparseMatrix<double> &PP,
                            const Eigen::SparseMatrix<double> &B,
                            const Eigen::SparseMatrix<double> &D,
                            const Eigen::SparseMatrix<double> &D_to_vel,
                            Eigen::VectorXd &u_new) {
    // solve linear system D * B * p = Δx*ρ/Δx * B * PP * u, i.e. ∇·∇p = ∇·u
    Eigen::SparseMatrix<double> A = B * D;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    std::cout << "compute success: " << (solver.info() == Eigen::Success) << "\n";
    Eigen::VectorXd p = solver.solve(dx * rho / dt * B * PP * u);
    std::cout << "solve success: " << (solver.info() == Eigen::Success) << "\n";
    Eigen::VectorXd dp = D_to_vel * D * p / dx;
    u_new = u - dt/rho * dp * p;
}