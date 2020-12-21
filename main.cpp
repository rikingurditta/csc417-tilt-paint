#include "include/paint_headers.h"
#include "make_cube.h"
#include <iostream>
#include <Eigen/Core>
#include <cstdlib>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/combine.h>
#include <random>

#include "pressure_projection_2d.h"
#include "fixed_zero_velocities_2d.h"
#include "div_matrix_2d.h"
#include "grad_matrix_2d.h"
#include "grad_grid_to_velocity_grid.h"
#include "advect_2d.h"
#include "apply_gravity_2d.h"
#include "particles_to_velocity_grid_2d.h"


bool rotate_grid(unsigned int key, Eigen::Matrix3d &R);

void apply_rotation(const Eigen::MatrixXd &R, const Eigen::RowVector3d &centre, const Eigen::MatrixXd &points,
                    Eigen::MatrixXd &result);

int main(int argc, char *argv[]) {
    // Main architecture goes as follows:
    // Values are hardcoded below, unless otherwise specified by command line arguments

    int nx = 100, ny = 100;
    double spacing = 0.02;
    double dt = 0.1;
    double rho = 1;
    // global gravity vector
    Eigen::Vector3d g(0, -9.8, 0);
    // current rotation of grid
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
    // original normal to 2d grid, don't need when/if we do 3d implementation
    Eigen::Vector3d n_original = Eigen::Vector3d(0, 0, 1);

    Eigen::Vector3d corner = Eigen::Vector3d(nx * spacing / -2.0, ny * spacing / -2.0, 0.0);
    Eigen::Vector3d centre = corner + Eigen::Vector3d(nx * spacing / 2, ny * spacing / 2, 0);

    Eigen::VectorXd num_particles_grid = Eigen::VectorXd::Zero(nx * ny);

    // create particles
    int n = 1000;
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> rand_x(0, nx);
    std::uniform_real_distribution<double> rand_y(0, ny);
    Eigen::MatrixXd particles = Eigen::MatrixXd::Zero(n, 2);
    Eigen::MatrixXd particle_velocities = Eigen::MatrixXd::Zero(n, 2);
    for (int i = 0; i < n; i++) {
        particles(i, 0) = rand_x(mt);
        particles(i, 1) = rand_y(mt);
        int px = floor(particles(i, 0));
        int py = floor(particles(i, 1));
        num_particles_grid(px + py * nx)++;
    }
    Eigen::Vector3d particles_centre = Eigen::Vector3d(nx * spacing / 2, ny * spacing / 2, 0);

    int vel_dx_grid_size = (nx + 1) * ny;
    int vel_dy_grid_size = nx * (ny + 1);
    Eigen::VectorXd u = Eigen::VectorXd::Zero(vel_dx_grid_size + vel_dy_grid_size);
    Eigen::VectorXd u_new = Eigen::VectorXd::Zero(vel_dx_grid_size + vel_dy_grid_size);
    // create matrix B which, given gradient, calculates divergence
    Eigen::SparseMatrix<double> B;
    div_matrix_2d(nx, ny, B);
    // create matrix PP which gets rid of boundary velocities
    Eigen::SparseMatrix<double> PP;
    fixed_zero_velocities(nx, ny, PP);
    // create gradient matrix D
    Eigen::SparseMatrix<double> D;
    grad_matrix_2d(nx, ny, D);

    Eigen::SparseMatrix<double> D_to_vel;
    grad_grid_to_velocity_grid(nx, ny, D_to_vel);

    std::cout << "B: " << B.rows() << ", " << B.cols() << "\n";
    std::cout << "PP: " << PP.rows() << ", " << PP.cols() << "\n";
    std::cout << "D: " << D.rows() << ", " << D.cols() << "\n";


    // main simulation loop
    auto simulate = [&](double delta_t) {
        advect_2d(delta_t, particles, particle_velocities);
        Eigen::Vector2d gravity;
        // current normal to grid
        Eigen::Vector3d n_curr = R * n_original;
        std::cout << "|n_curr|: " << n_curr.norm() << "\n";
        Eigen::Vector3d g_curr = R.transpose() * (g - g.dot(n_curr) * n_curr);
        Eigen::Vector2d g_curr_2d = g_curr.segment(0, 2);
        std::cout << "g_curr: " << g_curr.transpose() << "\n";
        apply_gravity_2d(delta_t, g_curr_2d, particle_velocities);
//        particles_to_velocity_grid_2d(particles, particle_velocities, nx, ny, u);
        pressure_projection_2d(u, nx, ny, spacing, delta_t, rho, PP, B, D, D_to_vel, u_new);
        u = u_new;
    };


    // visualizer code ----------------------------------------------------------------
    Eigen::RowVector3d center_viewer = Eigen::RowVector3d(0.0, 0.0, 0.0);
    // Helper function to rotate grid and update based on gravity

    // Define a grid by the origin, dimensions in x, y,
    // a corner, and the grid spacing
    Eigen::MatrixXd grid_points;
    make_grid(nx, ny, corner, spacing, grid_points);

    std::vector<Eigen::MatrixXd> Vs;
    std::vector<Eigen::MatrixXi> Fs;
    for (int i = 0; i < grid_points.rows(); i++) {
        if (num_particles_grid(i) > 0) {
            Eigen::MatrixXd V_cube;
            Eigen::MatrixXi F_cube;
            make_cube(grid_points.row(i), Eigen::Vector3d::Ones() * spacing, V_cube, F_cube);
            Vs.emplace_back(V_cube);
            Fs.emplace_back(F_cube);
        }
    }
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::combine(Vs, Fs, V, F);


    // Create a libigl Viewer object to toggle between point cloud and mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().show_lines = false;
    viewer.core().lighting_factor = 0;

    std::cout << R"(
    P,p      view point cloud
    pls replace commands here lol
    )";
    const auto set_points = [&]() {
        viewer.data().clear();
        // viewer.data().set_points(grid_points, center_viewer);
        viewer.data().set_mesh(V, F);
    };

    set_points();
    viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer &, unsigned int key, int) {
        switch (key) {
            case 'W':
            case 'w':
            case 'S':
            case 's':
            case 'A':
            case 'a':
            case 'D':
            case 'd':
            case 'R':
            case 'r':
                rotate_grid(key, R);
                break;
            case 'N':
            case 'n':
                simulate(dt);
        }
        Eigen::MatrixXd particles_3d = Eigen::MatrixXd::Zero(particles.rows(), 3);
        particles_3d.block(0, 0, particles.rows(), 2) = particles * spacing;
        Eigen::MatrixXd particles_3d_rot;
        apply_rotation(R, particles_centre, particles_3d, particles_3d_rot);
        viewer.data().set_points(particles_3d_rot, Eigen::RowVector3d::Ones());
        Eigen::MatrixXd V_rot;
        apply_rotation(R, centre, V, V_rot);
        viewer.data().set_vertices(V_rot);

        return false;
    };

    Eigen::MatrixXd C(F.rows(), 3);
    // TODO: rename colours
    Eigen::RowVector3d c1 = Eigen::RowVector3d(1.0, 1.0, 1.0);
    Eigen::RowVector3d c2 = Eigen::RowVector3d(1.0, 1.0, 0.0);
    Eigen::RowVector3d c3 = Eigen::RowVector3d(0.0, 1.0, 1.0);
    Eigen::RowVector3d c4 = Eigen::RowVector3d(1.0, 0.0, 1.0);
    Eigen::RowVector3d c5 = Eigen::RowVector3d(0.0, 0.0, 0.0);
    C << c1.replicate(F.rows() / 5, 1),
            c2.replicate(F.rows() / 5, 1),
            c3.replicate(F.rows() / 5, 1),
            c4.replicate(F.rows() / 5, 1),
            c5.replicate(F.rows() - 4 * (F.rows() / 5), 1);


    // paint_colours(nx, ny, C);

    viewer.data().set_colors(C);
    viewer.data().set_face_based(true);
    viewer.data().point_size = 2;
    viewer.launch();

    return EXIT_SUCCESS;
}

bool rotate_grid(unsigned int key, Eigen::Matrix3d &R) {
    float deg = 5.0 * M_PI / 180.0;
    Eigen::Matrix3d rot_right, rot_left, rot_up, rot_down;
    rot_right << 1.0, 0.0, 0.0,
            0.0, cos(deg), -sin(deg),
            0.0, sin(deg), cos(deg);
    rot_left << rot_right.transpose();

    rot_up << cos(deg), 0.0, -sin(deg),
            0.0, 1.0, 0.0,
            sin(deg), 0.0, cos(deg);
    rot_down << rot_up.transpose();

    switch (key) {
        case 'W':
        case 'w':
            std::cout << "up" << std::endl;
            R *= rot_up;
            return true;
        case 'S':
        case 's':
            std::cout << "down" << std::endl;
            R *= rot_down;
            return true;
        case 'A':
        case 'a':
            std::cout << "left" << std::endl;
            R *= rot_left;
            return true;
        case 'D':
        case 'd':
            std::cout << "right" << std::endl;
            R *= rot_right;
            return true;
        case 'R':
        case 'r':
            std::cout << "reset" << std::endl;
            R = Eigen::Matrix3d::Identity();
            return true;
    }
    return false;
}

void apply_rotation(const Eigen::MatrixXd &R, const Eigen::RowVector3d &centre, const Eigen::MatrixXd &points,
                    Eigen::MatrixXd &result) {
    Eigen::MatrixXd shift_centre = centre.replicate(points.rows(), 1);
    result = (points - shift_centre) * R.transpose() + shift_centre;
}