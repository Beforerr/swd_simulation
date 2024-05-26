# IDsHybridSimulation

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> IDsHybridSimulation

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "IDsHybridSimulation"
```
which auto-activate the project and enable local path handling from DrWatson.


## Installation (VPIC)

~~`vpic` needs `gcc` to run. `clang` and `oneapi` is not working though it successfully compiles.~~

## Installation (WarpX)

### Local Development (macOS)

- run `just env-macos` to setup the environment, then
- run `just py-warpx` to install `py-warpx`, for more details, see [justfile](justfile).

Notes:

- `OpenMPI requires both C and Fortran compilers!`
- `gcc` is not well supported in macOS: Installed by `spack` with `gcc` has bugs when compiling `openblas`.
    - [openblas 0.3.24 fails to build on aarch64-apple-darwin (all versions) in Homebrew · OpenMathLib/OpenBLAS](https://github.com/OpenMathLib/OpenBLAS/issues/4212) #issue
- warpx installed by `conda` has no `MPI` and `2d/3d` support.

```
spack env create warpx
spack env activate warpx
spack add py-warpx@develop ^warpx~openpmd
spack add ccache
spack concretize
for i in {1..4}; do nohup spack install >> install.txt 2>&1 & done
spack install
```

Using pip to build is slow

```
export BUILD_PARALLEL=8
export WARPX_MPI=ON
python3 -m pip wheel -v git+https://github.com/ECP-WarpX/WarpX.git
python3 -m pip install pywarpx*whl
```

### HPC

```
module swap hdf5 hdf5-mpi/1.12.2
module swap fftw fftw-mpi

spack add warpx%oneapi~openpmd dims="1,2,3"
```