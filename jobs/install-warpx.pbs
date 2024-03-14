#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=32:mpiprocs=1

source $HOME/.bashrc
micromamba activate warpx
micromamba list

module load cmake fftw-mpi
module swap hdf5 hdf5-mpi
module swap netcdf netcdf-mpi

cd $HOME/src/warpx
cmake -S . -B build \
  -DWarpX_DIMS="1;2;3" \
  -DWarpX_PYTHON=ON
cmake --build build --target pip_install -j 32