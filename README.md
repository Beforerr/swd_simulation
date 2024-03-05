

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