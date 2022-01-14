# maDGiCart-CH

![unit_test_ubuntu_20.04_github_runner_serial](https://github.com/mlohry/maDGiCart-CH/actions/workflows/unit_test_github_runner.yml/badge.svg?branch=master)
![unit_test_ubuntu_20.04_github_runner_openmp](https://github.com/mlohry/maDGiCart-CH/actions/workflows/unit_test_github_runner_openmp.yml/badge.svg?branch=master)
![unit_test_ubuntu_20.04_github_runner_cuda_](https://github.com/mlohry/maDGiCart-CH/actions/workflows/unit_test_github_runner_cuda.yml/badge.svg?branch=master)

## Running

After compilation, unit tests and a sample solution can be run with

    ./unit_testing
    ./maDGiCart

A longer running 3D example solution can be run with

    ./maDGiCart --dimension=3 --max_time_steps=10000 --eps2=0.0019974621629115655 --sigma=3.5306509073075123 --time_step_size=3e-5 --domain_x_begin=-3.14159 --domain_x_end=3.14159

If one chooses the `--time_integrator=ode23` option, error-adaptive time-stepping will be used. In addition, a global relative convergence tolerance can be specified with `--converged_rel_tol=value` whereby computation ends when absolute residual norm is less than this factor multiplied by the maximum residual over all steps. A complete run may then look like:

    ./maDGiCart --dimension=3 --max_time_steps=100000 --eps2=0.0019974621629115655 --sigma=3.5306509073075123 --domain_x_begin=-3.14159 --domain_x_end=3.14159 --time_integrator=ode23 --converged_rel_tol=1e-6

Alternatively one can use `--converged_abs_tol=value` e.g. `1e-5` for absolute convergence tolerance. 

## Plotting

Solutions are output to a VTK-standard `.vts` structured grid format. This can be loaded directly into programs such as Paraview, or with the included `plot_vts.py` script which utilizes the `pyvista` VTK frontend:

    python3 plot_vts.py outputfile.vts

Convergence history and configuration options are exported to `*.log` files which can be plotted with

    python3 /path/to/maDGiCart-CH/plot_history.py /path/to/logfile.log
     

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

### Singularity build from Docker

If you wish to build and run using Singularity, you can convert the Docker image (built as described above) using this workflow.

Get the image id from this command:

    docker images

Save that image as a tar file (say it has the Docker tag madg-cuda):

    docker save -t madg-cuda -o myimage.tar

Build the singularity image:

    singularity build myimage.sif docker-archive://myimage.tar

The resulting .sif image can be used with singularity. Note that you will need to set the environment variable LC_ALL=C. For example, you would start a singularity shell as follows:

    singularity shell --nv --env LC_ALL=C myimage.sif

Inside this shell, you would build using cmake/make in the same way as described in the Docker section.

## Benchmarking

System benchmarks can be run using `./benchmark_testing`. This measures the time spent in the evaluation of the right-hand-side of the discretized Cahn-Hilliard equation in 2D and 3D.

### Evaluations per second (higher is faster), 3D

AMD Ryzen 9 3900X 12-core CPU, nvidia GTX Titan Black GPU:

| Grid size     | Serial CPU  | OpenMP CPU, 12 threads | GTX Titan Black | 
| ------------- | ------------- | ------------- | ------------- |
| 64^3  | 84.48 | 385.71 | 353.88   |
| 128^3  | 10.15  | 41.52 |  48.66  |
| 256^3  | 1.20  | 4.59 |   5.78   |

AMD Ryzen Threadripper 1900X 8-Core CPU:

| Grid size     | Serial CPU  | OpenMP CPU, 8 threads |
| ------------- | ------------- | ------------- |
| 64^3  | 71.64 | 174.43| 
| 128^3  | 8.49  | 27.56 |  
| 256^3  | 0.93  | 2.94 |  


### Evaluations per second (higher is faster), 2D

AMD Ryzen 9 3900X 12-core CPU, nvidia GTX Titan Black GPU:

| Grid size     | Serial CPU  | OpenMP CPU, 12 threads | GTX Titan Black | 
| ------------- | ------------- | ------------- | ------------- |
| 64^2   | 11261.26  | 13812.15 | 8110.30  |
| 128^2  | 3045.99  | 7645.26 |  6684.49  |
| 256^2  | 744.82  | 2707.09 |   1194.89   |

AMD Ryzen Threadripper 1900X 8-Core CPU:

| Grid size     | Serial CPU  | OpenMP CPU, 8 threads |
| ------------- | ------------- | ------------- |
| 64^2  | 10582.01 | 11723.33 | 
| 128^2  | 2928.26  | 6540.22 |  
| 256^2  | 735.46 | 2123.14 |  
