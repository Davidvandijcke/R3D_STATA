# CLAUDE.md — R3D Stata Package

## What This Is

Stata port of the R3D R package implementing regression discontinuity design with distribution-valued outcomes. Two estimators: local polynomial (pointwise quantile-by-quantile) and Frechet regression (global). Multiplier bootstrap for uniform inference. Supports sharp and fuzzy designs.

**Repo:** `Davidvandijcke/R3D_STATA` on GitHub
**R package:** `Davidvandijcke/r3d` on GitHub
**Paper:** Van Dijcke (2025), "Regression Discontinuity Design with Distributional Outcomes"

## Build and Development

### Building the Plugin (Optional)
```bash
cd plugin/
make clean all
make install  # Copies plugin to ado/
```

### Running in Stata
```stata
# Run example
/Applications/Stata/StataSE.app/Contents/MacOS/stataSE -b do example.do

# Diagnose plugin
/Applications/Stata/StataSE.app/Contents/MacOS/stataSE -b do diagnose_plugin.do
```

### Installation
```stata
net install r3d, from("https://raw.githubusercontent.com/Davidvandijcke/R3D_STATA/main/") replace
```

## Architecture

Three-layer: ado (syntax) → Mata (algorithms) → Fortran plugin (optional LAPACK acceleration).

See `CODEBASE_MAP.md` for detailed structure.

## Key Algorithms

- **r3d_locpoly**: Local polynomial regression. Computes `(X'WX)^{-1} X'WY` for coefficients and `e_1' (X'WX)^{-1}` for intercept weights.
- **r3d_compute_quantiles**: Empirical quantiles using R type-7 interpolation.
- **r3d_estimate_density**: Silverman's rule `h = 1.06 * sd * n^{-1/5}` (matches R).
- **r3d_bandwidth_select**: MSE/IMSE-optimal bandwidth following paper Appendix A.3.
- **r3d_bootstrap**: Multiplier bootstrap for uniform confidence bands.

## Reference Implementation

The R package in `../` (same parent repo `r3d`) is the reference. All numerical results should match to machine precision where possible.

## Coding Guidelines

- `mata set matastrict on` — all variable declarations at function top
- Match R implementation numerically; when in doubt, follow R code
- Test against R-generated reference data in `tests/`
- Base plugin is optional; pure Mata must always work

## Git Workflow

```
feat/<description>    # New features
fix/<description>     # Bug fixes
docs/<description>    # Documentation
test/<description>    # Tests
```

Conventional commits. PR per issue.
