

## Installation (VPIC)

`vpic` needs `gcc` to run. `clang` and `oneapi` is not working though it successfully compiles. 

## Installation (WarpX)

### Local Development (macOS)

warpx installed by `conda` has no `MPI` support.

`OpenMPI requires both C and Fortran compilers!`

```
spack env create -d .
spack env activate .
spack add warpx~openpmd
spack add py-warpx
spack add openblas~fortran
spack concretize
for i in {1..4}; do nohup spack install >> install.txt 2>&1 & done
spack install
```

installed by `spack` with `gcc` has bugs when compiling `openblas`

- [openblas 0.3.24 fails to build on aarch64-apple-darwin (all versions) in Homebrew · OpenMathLib/OpenBLAS](https://github.com/OpenMathLib/OpenBLAS/issues/4212) #issue