#include "paint_headers.h"
#include <iostream>


// Make some paint colours

void paint_colours(
    const int nx,
    const int ny,
    Eigen::MatrixXd & colour_mat)
{
    // Colour matrices in a cool way here
    // diagonal
    int factor = 5;
    int mini_nx = nx / factor;
    int mini_ny = ny /factor;

    // I can't use block colours bc its 1d
    // unless i do it on another matrix then convert to 1d

    int o_x = 0;
    int o_y = 0;

    // one block
    for (int i = 0; i < mini_nx; i++ ) {
        for (int j = 0; i < mini_ny; j++) {
            int ind = nx * j + i;
            std::cout << ind <<"," << i <<"," << j<< std::endl;
            colour_mat.row(ind) = Eigen::RowVector3d(0.0, 0.1 * i, 0.2 * j);
        }
    }


    // for (int f = 0; f < factor; f++) {
    //     for (int i = o_x; i < o_x + factor; i++) {
    //         for (int j = o_y; i < o_y + factor; j++) {
    //             // Assign a colour 
    //             int ind = nx * j + i;
    //             std::cout << ind <<"," << o_x <<"," << o_y<< std::endl;
    //             // colour_mat.row(ind) = Eigen::RowVector3d(0.0, 0.1 * i, 0.2 * j);
    //         }
    //     }

    //     o_x += factor / 2;
    //     o_y += factor / 2;
    // }






}