# CSC417 Final Project - Tilt Paint

By Rikin Gurditta and Roya Shams-Zadeh-Amiri

---

## Video

Here is our video: https://drive.google.com/file/d/1ycEQbujZiCwDj8drVXzYqyF2lIv3_XMB/view?usp=sharing

## Report

Here is our report: [CSC417_Final_Project_Report.pdf](CSC417_Final_Project_Report.pdf)

## How to compile

Our code can be built like any other `cmake` project:

```bash
mkdir build
cd build
cmake ..
# or use below for release mode:
# cmake .. -DCMAKE_BUILD_TYPE=Release
./paint
```

## How to use

Upon running the program, if you cannot see any particles, you might need to pan upwards and rightwards to find them.

Use `W` and `S` to rotate the canvas around the x-axis and `A` and `D` to rotate it around the z-axis. Use `R` to reset the rotation.

Use `N` to step forward in time and perform one iteration of the simulation. usually the simulation will continue stepping if you hold `N` down.

## How to change parameters

Program parameters are all defined at the beginning of `main` if you would like to experiment with different parameters.