# maDGiCart-CH

![unit_test_ubuntu_20.04_github_runner](https://github.com/mlohry/maDGiCart-CH/actions/workflows/unit_test_github_runner.yml/badge.svg?branch=master)

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
    
In the cmake command line one can optionally specify one of `-DMADG_USE_SERIAL=On`, `-DMADG_USE_OPENMP=On`, or `-DMADG_USE_GPU=On` where serial is the default.

Then run unit tests and a sample solution with

    ./unit_testing
    ./maDGiCart


### Docker build - CPU

For consistent gcc CPU-based builds, a [Dockerfile](Dockerfile.gcc) is provided. The following workflow will build and run maDGiCart-CH in a container.

First, clone the repository on the host and make a build directory:

    git clone git@github.com:mlohry/maDGiCart-CH.git
    mkdir maDGiCart-CH-build

Build the image with tagged name `madg-gcc`:

    docker build -f ./maDGiCart-CH/Dockerfile.gcc  -t madg-gcc .

Start the image with the appropriate directories mounted and start a bash shell:

    docker run --rm -it --mount type=bind,source=$PWD/maDGiCart-CH,target=/maDGiCart-CH --mount type=bind,source=$PWD/maDGiCart-CH-build,target=/maDGiCart-CH-build madg-gcc bash

Compile as above:

    cd /maDGiCart-CH-build
    cmake /maDGiCart-CH -DMADG_USE_OPENMP=On  # optionally enable OpenMP
    make -j 8

Run unit tests and a sample solution:

    ./unit_testing
    ./maDGiCart
    