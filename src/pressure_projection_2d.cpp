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
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    solver.compute(A);
    std::cout << "compute success: " << (solver.info() == Eigen::Success) << "\n";
    Eigen::VectorXd p = solver.solve(dx * rho / dt * B * PP * u);
    std::cout << "solve success: " << (solver.info() == Eigen::Success) << "\n";
    Eigen::VectorXd dp = D_to_vel * D * p / dx;
    std::cout << "dp:\n";
    int vel_x_size = (nx + 1) * ny;
    int vel_y_size = nx * (ny + 1);
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx + 1; i++) {
            std::cout << dp(i + j * (nx + 1)) << " ";
        }
        std::cout << "\n";
    }
    for (int j = 0; j < ny + 1; j++) {
        for (int i = 0; i < nx; i++) {
            std::cout << dp(vel_x_size + i + j * nx) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n\n";
    u_new = u - dt/rho * dp;
}