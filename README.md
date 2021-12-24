# maDGiCart-CH

![unit_test_ubuntu_20.04_github_runner](https://github.com/mlohry/maDGiCart-CH/actions/workflows/unit_test_github_runner.yml/badge.svg?branch=master)

## Building

Using an out-of-source build, use cmake:

    cmake ./path/to/maDiCart-CH
    make -j
    
In the cmake command line one can optionally specify one of `-DMADG_USE_SERIAL=On`, `-DMADG_USE_OPENMP=On`, or `-DMADG_USE_GPU=On` where serial is the default.

Then run unit tests and a sample solution with

    ./unit_testing
    ./maDGiCart
