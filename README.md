

## Installation (VPIC)

`vpic` needs `gcc` to run. `clang` and `oneapi` is not working though it successfully compiles. 

## Installation (WarpX)

### Local Development (macOS)

warpx installed by `conda` has no `MPI` and `2d/3d` support.

Notes:

- `OpenMPI requires both C and Fortran compilers!`
- `gcc` is not well supported in macOS: Installed by `spack` with `gcc` has bugs when compiling `openblas`.
    - [openblas 0.3.24 fails to build on aarch64-apple-darwin (all versions) in Homebrew · OpenMathLib/OpenBLAS](https://github.com/OpenMathLib/OpenBLAS/issues/4212) #issue

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

## HPC

```
module swap hdf5 hdf5-mpi/1.12.2
module swap fftw fftw-mpi

spack add warpx%oneapi~openpmd dims="1,2,3"
```