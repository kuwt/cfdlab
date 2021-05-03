This repository is a CFD Solver developed by Qianqian, Christina and Wing To from TUM CSE. It is built on top of Fluidchen from TUM Informatics, Chair of Scientific Computing in Computer Science

## Software Requirements
* VTK 7 or higher
* GCC 9 (optional)

### Setup of VTK and GCC 9 (Ubuntu **20.04**)

```shell
apt-get update &&
apt-get upgrade -y &&
apt-get install -y build-essential cmake libvtk7-dev libfmt-dev
```

## Installing

```shell
git clone https://gitlab.lrz.de/cfdlab_2021/fluidchen-skeleton.git GroupC_CFDLab
cd GroupC_CFDLab
mkdir build && cd build
cmake ..
make
```

If you are running GCC 8 or lower, you have to set gpp9 to false in CMakeList.txt. If you have many gcc installed into your system. It is better to use -DCMAKE_C_COMPILER option to specify the gcc you want to use. For example, to use gcc8, replace the cmake .. line above to :
```shell
cmake -DCMAKE_C_COMPILER=/usr/bin/gcc-8 -DCMAKE_CXX_COMPILER=/usr/bin/g++-8  ..
```

If you want to debug using gdb, do not forget to add debug option.
```shell
cmake -DCMAKE_BUILD_TYPE=Debug ..
```

## Running

In order to run the code, the case file should be given as input parameter. The lid driven cavity case files are located in the `example_cases` directory. In the build directory, 

```shell
./fluidchen /full_path_to_this_repo_in_your_system/example_cases/LidDrivenCavity/LidDrivenCavity.dat
```

This will run the case file and create the output folder `/full_path_to_this_repo_in_your_system/example_cases/LidDrivenCavity/LidDrivenCavityOutput` which holds the `.vtk` files of the solution.  Note that this may require write permissions in the given directory.

## Result Visualization

You can use paraview to visualize the output .vtk file. Here is an example showing the velocity magnitude of the LidDrivenCavity problem.

![sample LidDrivenCavity Velocity](docs/images/LidDrivenCavitySampleVelocityImage.png "Velocity")
