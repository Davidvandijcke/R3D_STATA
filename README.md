# R3D: Regression Discontinuity with Distribution-valued Outcomes

Stata implementation of the R3D estimator from Van Dijcke (2025), "Regression Discontinuity Design with Distributional Outcomes."

Estimates treatment effects on entire outcome distributions at an RD cutoff, rather than just the mean. Supports sharp and fuzzy designs, two estimation methods, multiplier bootstrap for uniform inference, and data-driven bandwidth selection.

Based on the [R package](https://github.com/Davidvandijcke/r3d).

## Installation

### From GitHub (recommended)

```stata
net install r3d, from("https://raw.githubusercontent.com/Davidvandijcke/R3D_STATA/main/") replace
```

### Manual installation

1. Download or clone this repository
2. Copy `ado/*.ado`, `mata/*.mata`, and `doc/*.sthlp` to your Stata PERSONAL or PLUS directory
3. Optionally copy the compiled plugin (`ado/r3d_plugin.plugin`) for faster computation

Find your PERSONAL directory with:
```stata
sysdir
```

## Quick Start

```stata
// Generate example data: 500 obs, 3 outcome draws per unit
set obs 500
set seed 12345
gen x = rnormal()
forvalues i = 1/3 {
    gen y`i' = 1 + 0.5*(x >= 0) + 0.2*x + rnormal(0, 1)
}

// Sharp RD at cutoff = 0
r3d x y1 y2 y3, cutoff(0)

// With bandwidth selection and bootstrap inference
r3d x y*, cutoff(0) bwselect bootstrap(999)

// Fuzzy RD
gen treat = (x >= 0) * (runiform() > 0.2)
r3d x y*, cutoff(0) fuzzy(treat)
```

## Syntax

```
r3d runvar outcomevar1 [outcomevar2 ...] [if] [in],
    cutoff(#) [method(simple|frechet)] [fuzzy(varname)]
    [polynomial(#)] [bandwidth(numlist)] [bwselect]
    [kernel(triangular|epanechnikov|uniform)]
    [quantiles(numlist)] [nquantiles(#)]
    [bootstrap(#)] [level(#)]
    [weights(varname)] [nograph] [saving(filename)] [replace]
```

### Key options

| Option | Description | Default |
|--------|-------------|---------|
| `cutoff(#)` | RD cutoff value | *required* |
| `method()` | `simple` (pointwise) or `frechet` (global) | `simple` |
| `fuzzy()` | Treatment variable for fuzzy RD | sharp RD |
| `polynomial(#)` | Local polynomial order | 2 |
| `bandwidth()` | Bandwidth(s) for local regression | data-driven |
| `bwselect` | Automatic bandwidth selection | off |
| `kernel()` | Kernel function | `triangular` |
| `nquantiles(#)` | Number of quantile grid points | 99 |
| `bootstrap(#)` | Number of bootstrap replications | no bootstrap |

## Methods

**Simple** (`method(simple)`): Estimates treatment effects pointwise at each quantile using separate local polynomial regressions. Fast and flexible.

**Frechet** (`method(frechet)`): Uses global Frechet regression in Wasserstein space with isotonic projection to ensure monotonicity of estimated quantile functions. More efficient when the monotonicity constraint binds.

## Stored results

After estimation, `r3d` stores:

| Result | Description |
|--------|-------------|
| `e(b)` | Treatment effect vector (by quantile) |
| `e(V)` | Variance-covariance matrix |
| `e(quantiles)` | Quantile grid |
| `e(N)` | Number of observations |
| `e(cutoff)` | RD cutoff |
| `e(method)` | Estimation method |
| `e(h_num)` | Bandwidth for outcome |
| `e(h_den)` | Bandwidth for treatment (fuzzy only) |

## Architecture

The package has three layers:

1. **Ado files** (`ado/`): Stata command syntax, option parsing, result storage
2. **Mata** (`mata/`): Core algorithms — quantile computation, local polynomial regression, bandwidth selection, bootstrap, isotonic regression
3. **Fortran plugin** (`plugin/`, optional): LAPACK-accelerated local polynomial solver for performance

The Fortran plugin is optional. If not available, pure Mata is used automatically.

### Building the plugin

Requires `gcc`, `gfortran`, and LAPACK. See `plugin/README.md` for details.

```bash
cd plugin/
make clean all
make install   # copies to ado/
```

## Post-estimation commands

```stata
// Bootstrap inference (after r3d)
r3d_bootstrap, reps(999) tests(nullity homogeneity)

// Standalone bandwidth selection
r3d_bwselect x y*, cutoff(0)
```

## References

Van Dijcke, D. (2025). "Regression Discontinuity Design with Distributional Outcomes." *Working Paper*.

Cattaneo, M. D., Idrobo, N., and Titiunik, R. (2019). *A Practical Introduction to Regression Discontinuity Designs: Foundations*. Cambridge University Press.

## Author

David Van Dijcke, University of Michigan
dvdijcke@umich.edu

## License

MIT
