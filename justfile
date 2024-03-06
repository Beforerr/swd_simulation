home_dir := env_var('HOME')
vpic := home_dir + "/src/vpic-kokkos/build/bin/vpic"

compile:
  cd {{invocation_directory()}}; {{vpic}} *.cxx

q-info:
  qhist -u $USER
  qstat -u $USER

env-install:
   micromamba env create --file environment.yml

env-update:
   micromamba install --file environment.yml

spack-compilers:
   code $HOME/.spack/darwin/compilers.yaml

warpx:
   module load cmake
   cd $HOME/src/warpx
   cmake -S . -B build \
      -DWarpX_DIMS="1;2;3" \
      -DPython_EXECUTABLE=$(which python3) \
      -DWarpX_PYTHON=ON 
   cmake --build build -j 8
   cmake --build build --target pip_install -j 8
   # python3 -m pip wheel -v .
   # python3 -m pip install pywarpx*whl

install-warpx-ncar:
   #!/bin/bash
   # DEBUG: not working
   spack install warpx%dpcpp


# not working
copy:
   ssh-copy-id zijin@derecho.hpc.ucar.edu