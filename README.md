

## Installation

### Local Development (macOS)

warpx installed by `conda` has no `MPI` support.

```
spack env create warpx
spack env activate warpx
spack add warpx%gcc
spack install
# despacktivate && spack env remove warpx-gcc
```

installed by `spack` with `gcc` has bugs when compiling `openblas`

- [openblas 0.3.24 fails to build on aarch64-apple-darwin (all versions) in Homebrew · OpenMathLib/OpenBLAS](https://github.com/OpenMathLib/OpenBLAS/issues/4212) #issue