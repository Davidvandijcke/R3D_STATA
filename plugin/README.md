# R3D Plugin Build Instructions

Optional Fortran/LAPACK plugin for accelerated local polynomial regression.

## Prerequisites

- **gcc** (C compiler)
- **gfortran** (Fortran compiler)
- **LAPACK/BLAS** (macOS: Accelerate framework; Linux: `liblapack-dev libblas-dev`)
- **Stata plugin SDK** (`stplugin.h`, `stplugin.c` — included in this directory)

## Build

```bash
cd plugin/
make clean all
make install   # copies plugin to ../ado/
```

## Platform Notes

**macOS (Apple Silicon / Intel):**
- Uses `-framework Accelerate` for LAPACK
- Static-links gfortran/quadmath to avoid runtime dependencies
- Tested with gfortran 15.x (Homebrew)

**Linux:**
- Uses `-llapack -lblas -lgfortran`
- Install: `sudo apt-get install gfortran liblapack-dev libblas-dev`

**Windows:**
- Uses `-llapack -lblas -lgfortran`
- Requires MinGW or similar GCC toolchain

## Architecture

```
r3d_plugin.c          C wrapper implementing Stata plugin API
  ↓ calls
r3d_fortran.f90       Fortran module with locweights() subroutine
  ↓ uses
LAPACK (DGETRF/DGETRS)  LU factorization and solve
```

The plugin computes `(X'WX)^{-1}X'WY` and `e_1'(X'WX)^{-1} * basis * K`
for each quantile, using LAPACK for numerical stability.

## Limits

- Maximum 500 quantiles (compile-time limit in `r3d_plugin.c`)
- Polynomial order 0-9
