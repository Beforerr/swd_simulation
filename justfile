# not working
env-install:
   micromamba env create --file environment.yml

env-update:
   micromamba install --file environment.yml

copy:
   ssh-copy-id zijin@derecho.hpc.ucar.edu

spack-compilers:
   code $HOME/.spack/darwin/compilers.yaml