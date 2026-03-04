# R3D Stata — Codebase Map

## Directory Structure

```
R3D_STATA/
├── ado/                          # Stata command files
│   ├── r3d.ado                   # Main estimation command (syntax parsing, orchestration)
│   ├── r3d_bootstrap.ado         # Post-estimation bootstrap command
│   ├── r3d_bwselect.ado          # Standalone bandwidth selection command
│   ├── r3d_plugin_load.ado       # Platform-specific plugin loader with Mata fallback
│   └── r3d_setup.ado             # Post-install setup helper
├── mata/                         # Mata computational layer
│   ├── r3d_mata.mata             # Core algorithms (~1555 lines):
│   │                               - r3d_kernel(): kernel function implementations
│   │                               - r3d_variance(): sample variance
│   │                               - r3d_compute_quantiles(): empirical quantile computation
│   │                               - r3d_estimate_density(): density at cutoff
│   │                               - r3d_kernel_matrices(): bias/variance kernel matrices
│   │                               - r3d_fit_global_poly(): global polynomial pilot
│   │                               - r3d_bandwidth_select(): MSE/IMSE-optimal bandwidth
│   │                               - r3d_prepare_bandwidth_matrix(): bandwidth validation
│   │                               - r3d_locpoly(): local polynomial regression
│   │                               - r3d_isotonic(): PAVA isotonic regression
│   │                               - r3d_estimate_core(): main estimation engine
│   │                               - r3d_simple(): simple method wrapper
│   │                               - r3d_frechet(): Frechet method wrapper
│   │                               - r3d_bootstrap(): multiplier bootstrap
│   │                               - r3d_gini_from_quantile(): Gini from quantile function
│   │                               - r3d_test_gini(): Gini difference test
│   └── r3d_plugin_interface.mata # Bridge between Mata and C/Fortran plugin:
│                                   - r3d_plugin_available(): check plugin status
│                                   - r3d_plugin_locweights(): call Fortran via Stata vars
├── plugin/                       # Optional C/Fortran plugin
│   ├── r3d_plugin.c              # C wrapper for Stata plugin API
│   ├── r3d_fortran.f90           # Fortran locweights with LAPACK (DGETRF/DGETRS)
│   ├── Makefile                  # Cross-platform build (macOS/Linux/Windows)
│   ├── stplugin.c                # Stata plugin SDK source
│   └── stplugin.h                # Stata plugin SDK header
├── doc/                          # Help files
│   ├── r3d.sthlp                 # Main command help
│   ├── r3d_bootstrap.sthlp       # Bootstrap help
│   └── r3d_bwselect.sthlp        # Bandwidth selection help
├── tests/                        # Test suite
│   └── (to be populated)
├── example.do                    # End-to-end usage example
├── diagnose_plugin.do            # Plugin loading diagnostic
├── install_r3d.do                # Local installation script
├── verify_package.do             # Package structure verification
├── stata.toc                     # Package index
├── r3d.pkg                       # Package manifest
├── r3d.ancillary                 # Ancillary file list
├── LICENSE                       # MIT License
└── CLAUDE.md                     # AI assistant context
```

## Architecture

Three-layer design:
1. **Stata ado** — User-facing commands, syntax parsing, result display
2. **Mata** — Core numerical algorithms (auto-compiled on first use)
3. **Fortran plugin** — Optional accelerated local polynomial (LAPACK-based)

The plugin is optional; the package falls back to pure Mata if unavailable.

## Data Flow

```
r3d.ado (parse syntax, prep data)
  → r3d_compute_quantiles() (empirical quantiles from Y vars)
  → r3d_bandwidth_select() (MSE-optimal bandwidths)
  → r3d_simple() / r3d_frechet()
      → r3d_locpoly() or r3d_plugin_locweights()
      → r3d_isotonic() (Frechet only)
  → r3d_bootstrap() (optional, multiplier bootstrap)
  → ereturn post (store results)
```
