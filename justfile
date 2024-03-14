home_dir := env_var('HOME')
vpic := home_dir + "/src/vpic-kokkos/build/bin/vpic"
project := "UCLA0040"
queue := "main"

compile:
  cd {{invocation_directory()}}; {{vpic}} *.cxx

q-sub:
  cd {{invocation_directory()}}; qsub -A {{project}} -q {{queue}}  *.pbs

q-subf file:
  cd {{invocation_directory()}}; qsub -A {{project}} -q {{queue}}  {{file}}

q-info:
  qhist -u $USER
  -qstat -u $USER

env-install file="environment.yml":
  micromamba env create --file {{file}}

env-update file="environment.yml":
  micromamba install --file {{file}}

spack-compilers:
  code $HOME/.spack/darwin/compilers.yaml

warpx:
  cd $HOME/src/WarpX
  cmake -S . -B build
  cmake --build build -j 8

py-warpx:
  #!/bin/sh
  cd $HOME/src/WarpX
  cmake -S . -B build \
    -DWarpX_DIMS="1;2;3" \
    -DWarpX_PYTHON=ON \
    -DWarpX_MPI=OFF
  cmake --build build --target pip_install -j 8

install-warpx-ncar:
  #!/bin/bash
  # DEBUG: not working
  spack install warpx%dpcpp

# not working
copy:
  ssh-copy-id zijin@derecho.hpc.ucar.edu

clean:
  #!/usr/bin/env bash
  cd {{invocation_directory()}}
  rm -rf *\.{e,o}*
  rm -f *.dpkl
  rm -f Backtrace.*
  rm -f warpx_used_inputs
  rm -rf diags