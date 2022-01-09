# maDGiCart-CH

![unit_test_ubuntu_20.04_github_runner_serial](https://github.com/mlohry/maDGiCart-CH/actions/workflows/unit_test_github_runner.yml/badge.svg?branch=master)
![unit_test_ubuntu_20.04_github_runner_openmp](https://github.com/mlohry/maDGiCart-CH/actions/workflows/unit_test_github_runner_openmp.yml/badge.svg?branch=master)
![unit_test_ubuntu_20.04_github_runner_cuda_](https://github.com/mlohry/maDGiCart-CH/actions/workflows/unit_test_github_runner_cuda.yml/badge.svg?branch=master)

## Compiling

### CMake Options

The primary cmake command line options are

    -DMADG_USE_SERIAL=On
    -DMADG_USE_OPENMP=On
    -DMADG_USE_GPU=On

where `-DMADG_USE_SERIAL=On` is the default. Only one of these options can be `On` for a build.


### Ubuntu 20.04

On Ubuntu 20.04 the pre-requisites for building can be installed by

    sudo apt install libboost-log-dev libboost-program-options-dev libboost-regex-dev libboost-thread-dev libboost-filesystem-dev cmake g++

where other OS's should be similar. Using an out-of-source build, use cmake:

    cmake ./path/to/maDiCart-CH
    make -j

Then run unit tests and a sample solution with

    ./unit_testing
    ./maDGiCart

A longer running 3D example solution can be run with

    ./maDGiCart --dimension=3 --max_time_steps=10000 --eps2=0.0019974621629115655 --sigma=3.5306509073075123 --time_step_size=3e-5 --domain_x_begin=-3.14159 --domain_x_end=3.14159

### Docker build

For consistent builds, a [Dockerfile for CPU gcc builds](Dockerfile.gcc) and a [Dockerfile for GPU CUDA builds](Dockerfile.cuda) are provided. The following workflow will build and run maDGiCart-CH in a container.

One time only (assuming the docker images are cached), build either or both docker images for the CPU or GPU builds:

    docker build -f ./maDGiCart-CH/Dockerfile.gcc  -t madg-gcc .
    docker build -f ./maDGiCart-CH/Dockerfile.cuda  -t madg-cuda .

Clone the repository on the host and make a build directory:

    git clone git@github.com:mlohry/maDGiCart-CH.git
    mkdir maDGiCart-CH-build

Start the image with the appropriate directories mounted and start a bash shell for the CPU build:

    docker run --rm -it \
    --mount type=bind,source=$PWD/maDGiCart-CH,target=/maDGiCart-CH \
    --mount type=bind,source=$PWD/maDGiCart-CH-build,target=/maDGiCart-CH-build \
    madg-gcc bash
    
or for the GPU build (note the --gpus all option)

    docker run --gpus all --rm -it \
    --mount type=bind,source=$PWD/maDGiCart-CH,target=/maDGiCart-CH \
    --mount type=bind,source=$PWD/maDGiCart-CH-build,target=/maDGiCart-CH-build \
    madg-cuda bash

Compile as above for CPU:

    cd /maDGiCart-CH-build
    cmake /maDGiCart-CH -DMADG_USE_OPENMP=On # or -DMADG_USE_SERIAL=On
    make -j 8
    
or for GPU:

    cd /maDGiCart-CH-build
    cmake /maDGiCart-CH -DMADG_USE_GPU=On
    make -j 8

Run unit tests and a sample solution:

    ./unit_testing
    ./maDGiCart


## Plotting

Solutions are output to a VTK-standard `.vts` structured grid format. This can be loaded directly into programs such as Paraview, or with the included `plot_vts.py` script which utilizes the `pyvista` VTK frontend:

    python3 plot_vts.py outputfile.vts

