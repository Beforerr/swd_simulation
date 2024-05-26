import 'files/quarto.just'
overleaf_repo := ""
home_dir := env_var('HOME')

default:
  just --list

# install/update warpx
install-warpx:
  just py-warpx

warpx:
  cd $HOME/src/WarpX
  cmake -S . -B build
  cmake --build build -j 8

py-warpx:
  #!/bin/sh
  mkdir -p $HOME/src/WarpX
  cd $HOME/src/WarpX
  git pull # update the source code
  cmake -S . -B build \
    -DWarpX_DIMS="1;2;3" \
    -DWarpX_PYTHON=ON
  cmake --build build --target pip_install -j 8

clean:
  #!/usr/bin/env bash
  cd {{invocation_directory()}}
  rm -rf *\.{e,o}*
  rm -f *.dpkl
  rm -f Backtrace.*
  rm -rf diags
  find . -type f -name 'warpx_used_inputs' -exec rm {} +

run:
  #!/usr/bin/env bash
  cd {{invocation_directory()}}
  ipython inputs.ipynb

picviewer:
  cd dim_2_beta_0.25_theta_60/diags && picviewer
  cd dim_3_beta_0.25_theta_60/diags && picviewer