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

env-macos:
  micromamba env create --file environment.deps.yml
  just env-update

spack-compilers:
  code $HOME/.spack/darwin/compilers.yaml

warpx:
  cd $HOME/src/WarpX
  cmake -S . -B build
  cmake --build build -j 8

py-warpx:
  #!/bin/sh
  cd $HOME/src/WarpX
  git pull # update the source code
  cmake -S . -B build \
    -DWarpX_DIMS="1;2;3" \
    -DWarpX_PYTHON=ON
  cmake --build build --target pip_install -j 8

install-warpx-ncar:
  #!/bin/bash
  # DEBUG: not working
  spack install warpx%dpcpp

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

preview:
  quarto preview --no-render

publish:
  quarto publish gh-pages --no-render --no-prompt

vpic:
  micromamba env create vtk pyvista pyqt VisualPIC --name vpic
  micromamba run -n vpic vpic -h
  pipx install VisualPIC --preinstall pyqt5 --preinstall vtk --preinstall pyvista
  pipx inject visualpic pyqt5 vtk pyvista
  vpic dim_2_beta_0.25_theta_60/diags/diag1/ -Bx -By -Bz
  vpic dim_2_beta_0.25_theta_60/diags/diag1/ -Jx -Jy -Jz

  micromamba run -n vpic vpic3d dim_3_beta_0.25_theta_60/diags/diag1/ -Bx -By -Bz
  micromamba run -n vpic vpic3d dim_3_beta_0.25_theta_60/diags/diag1/ -Jx
  micromamba run -n vpic vpic3d dim_3_beta_0.25_theta_60/diags/diag1/ -Jz


picviewer:
  cd dim_2_beta_0.25_theta_60/diags && picviewer
  cd dim_3_beta_0.25_theta_60/diags && picviewer