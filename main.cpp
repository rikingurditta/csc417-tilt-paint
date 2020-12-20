#include "include/paint_headers.h"
#include "make_cube.h"
#include <iostream>
#include <Eigen/Core>
#include <cstdlib>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/combine.h>


bool rotate_grid(int key, Eigen::Matrix3d & R) {
    float deg = 5.0 * M_PI / 180.0;
    Eigen::Matrix3d rot_right, rot_left, rot_up, rot_down;
    rot_right << 1.0, 0.0, 0.0,
                0.0, cos(deg), -sin(deg),
                0.0, sin(deg), cos(deg); 
    rot_left << 1.0, 0.0, 0.0,
                0.0, cos(-deg), -sin(-deg),
                0.0, sin(-deg), cos(-deg); 

    rot_up << cos(deg), 0.0, -sin(deg),
                0.0, 1.0, 0.0,
                sin(deg), 0.0, cos(deg);
    rot_down << cos(-deg), 0.0, -sin(-deg),
                0.0, 1.0, 0.0,
                sin(-deg), 0.0, cos(-deg);

    switch(key)
        {
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
        }
    return false;
}



int main(int argc, char *argv[]) {
    // Main architecture goes as follows:
    // Values are hardcoded below, unless otherwise specified by command line arguments

    int nx = 100, ny = 100;
    Eigen::RowVector3d center_viewer = Eigen::RowVector3d(0.0, 0.0, 0.0);
    double spacing = 0.02;
    Eigen::Vector3d corner = Eigen::Vector3d(nx * spacing/-2.0, ny * spacing/-2.0, 0.0);

    // Helper function to rotate grid and update based on gravity
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();




    // Define a grid by the origin, dimensions in x, y, 
    // a corner, and the grid spacing
    Eigen::MatrixXd grid_points;
    make_grid(nx, ny, corner, spacing, grid_points);

    std::vector<Eigen::MatrixXd> Vs;
    std::vector<Eigen::MatrixXi> Fs;
    for (int i = 0; i < grid_points.rows(); i++) {
        Eigen::MatrixXd V_cube;
        Eigen::MatrixXi F_cube;
        make_cube(grid_points.row(i), Eigen::Vector3d::Ones() * spacing, V_cube, F_cube);
        Vs.emplace_back(V_cube);
        Fs.emplace_back(F_cube);
    }
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::combine(Vs, Fs, V, F);


    // Create a libigl Viewer object to toggle between point cloud and mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().show_lines = false;

    std::cout<<R"(
    P,p      view point cloud
    pls replace commands here lol
    )";
    const auto set_points = [&]()
    {
        viewer.data().clear();
        // viewer.data().set_points(grid_points, center_viewer);
        viewer.data().set_mesh(V, F);


    };
    set_points();
    viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer&, unsigned int key,int)
    {
        rotate_grid(key, R);
        // viewer.data().set_points(V * R.transpose(), center_viewer);
        viewer.data().set_vertices(V * R.transpose());

        return false;
    };
    viewer.data().point_size = 2;
    viewer.launch();
    
    std::cout << "Hello, World!" << std::endl;
    return EXIT_SUCCESS;
}
