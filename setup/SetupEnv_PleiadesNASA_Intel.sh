#!/bin/bash

CWD=`pwd`
export COMPILER_ID=Intel
export FC=mpif90
export CC="icc -lmpi"
export CXX="icpc -lmpi"
export FFTW_PATH=/home5/aghate/PadeOps/dependencies/fftw-3.3.5
export DECOMP_PATH=/home5/aghate/PadeOps/dependencies/2decomp_fft
export VTK_IO_PATH=/home5/aghate/PadeOps/dependencies/Lib_VTK_IO/build
export HDF5_PATH=/home5/aghate/PadeOps/dependencies/hdf5-1.8.18
