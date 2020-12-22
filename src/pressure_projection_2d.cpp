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
    // solve linear system D * B * p = Δx*ρ/Δt * B * PP * u, i.e. Δt/ρ ∇·∇p = ∇·u
    Eigen::SparseMatrix<double> A = B * D;
    // using QR solver instead of a conjugate gradient solver (such as BiCGSTAB) because the CG solvers weren't working
    // this is probably a hint that we are doing something wrong but it is too hard to figure out how to fix it
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    solver.compute(A);
    // std::cout << "compute success: " << (solver.info() == Eigen::Success) << "\n";
    Eigen::VectorXd p = solver.solve(dx * rho / dt * B * PP * u);
    // std::cout << "solve success: " << (solver.info() == Eigen::Success) << "\n";
    Eigen::VectorXd dp = D_to_vel * D * p / dx;
    // update velocity to make it divergence free
    u_new = u - dt/rho * dp;
}