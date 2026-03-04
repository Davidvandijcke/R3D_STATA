/*
 * test_equivalence_core.do
 * Tests numerical equivalence between Stata R3D and R R3D package.
 *
 * Prerequisites:
 *   1. Run `Rscript generate_reference_data.R` first to create ref_*.csv files
 *   2. R3D Stata package installed or Mata compiled
 *
 * Usage:
 *   cd tests/
 *   do test_equivalence_core.do
 */

clear all
set more off

// Compile Mata if not already compiled
quietly capture mata: mata which r3d_compute_quantiles
if _rc != 0 {
    quietly do "../mata/r3d_mata.mata"
}
quietly capture mata: mata which r3d_plugin_available
if _rc != 0 {
    quietly do "../mata/r3d_plugin_interface.mata"
}

local tolerance_coef = 1e-6
local tolerance_bw   = 1e-4
local n_pass = 0
local n_fail = 0
local n_total = 0

program define report_test
    args test_name result max_diff tolerance
    local n_total = ${n_total} + 1
    global n_total = `n_total'
    if `result' == 1 {
        di as text "  PASS: `test_name' (max diff = " as result %12.2e `max_diff' ///
           as text " < " as result %8.2e `tolerance' as text ")"
        local n_pass = ${n_pass} + 1
        global n_pass = `n_pass'
    }
    else {
        di as error "  FAIL: `test_name' (max diff = " as result %12.2e `max_diff' ///
           as text " >= " as result %8.2e `tolerance' as text ")"
        local n_fail = ${n_fail} + 1
        global n_fail = `n_fail'
    }
end

// ============================================================================
// Case 1: Sharp-Simple, nq=20, Epanechnikov
// ============================================================================
di _n as text "{hline 70}"
di as text "Case 1: Sharp-Simple, nq=20, Epanechnikov"
di as text "{hline 70}"

// Load data
import delimited "ref_data_sharp_simple_20.csv", clear

// Get variable list for Q columns
unab qvars : Q*
local nq : word count `qvars'

// Generate quantile grid
local quantiles
forvalues i = 1/`nq' {
    local q = `i' / (`nq' + 1)
    local quantiles `quantiles' `q'
}

// Run R3D
r3d X `qvars', cutoff(0) method(simple) polynomial(2) kernel(epanechnikov) ///
    quantiles(`quantiles') nograph

// Load reference results
preserve
import delimited "ref_results_sharp_simple_20.csv", clear

// Compare tau (treatment effects)
local max_diff_tau = 0
forvalues i = 1/`nq' {
    local ref_tau = tau[`i']
    local stata_tau = _b[q`i']
    local diff = abs(`stata_tau' - `ref_tau')
    if `diff' > `max_diff_tau' local max_diff_tau = `diff'
}
local pass = (`max_diff_tau' < `tolerance_coef')
report_test "Sharp-Simple tau" `pass' `max_diff_tau' `tolerance_coef'

// Compare bandwidths
local max_diff_bw = 0
forvalues i = 1/`nq' {
    local ref_bw = bw_num[`i']
    mata: st_numscalar("__bw", st_matrix("e(bandwidth_num)")[`i',1])
    local stata_bw = scalar(__bw)
    local diff = abs(`stata_bw' - `ref_bw')
    if `diff' > `max_diff_bw' local max_diff_bw = `diff'
}
local pass = (`max_diff_bw' < `tolerance_bw')
report_test "Sharp-Simple bandwidths" `pass' `max_diff_bw' `tolerance_bw'

// Compare intercepts
local max_diff_int = 0
forvalues i = 1/`nq' {
    local ref_intp = int_plus[`i']
    mata: st_numscalar("__ip", st_matrix("e(int_plus)")[1,`i'])
    local stata_intp = scalar(__ip)
    local diff = abs(`stata_intp' - `ref_intp')
    if `diff' > `max_diff_int' local max_diff_int = `diff'
}
local pass = (`max_diff_int' < `tolerance_coef')
report_test "Sharp-Simple int_plus" `pass' `max_diff_int' `tolerance_coef'
restore

// ============================================================================
// Case 2: Sharp-Simple, nq=99, Triangular
// ============================================================================
di _n as text "{hline 70}"
di as text "Case 2: Sharp-Simple, nq=99, Triangular"
di as text "{hline 70}"

import delimited "ref_data_sharp_simple_99.csv", clear
unab qvars : Q*
local nq : word count `qvars'

local quantiles
forvalues i = 1/`nq' {
    local q = `i' / 100
    local quantiles `quantiles' `q'
}

r3d X `qvars', cutoff(0) method(simple) polynomial(2) kernel(triangular) ///
    quantiles(`quantiles') nograph

preserve
import delimited "ref_results_sharp_simple_99.csv", clear
local max_diff_tau = 0
forvalues i = 1/`nq' {
    local ref_tau = tau[`i']
    local stata_tau = _b[q`i']
    local diff = abs(`stata_tau' - `ref_tau')
    if `diff' > `max_diff_tau' local max_diff_tau = `diff'
}
local pass = (`max_diff_tau' < `tolerance_coef')
report_test "Sharp-Simple-99q tau" `pass' `max_diff_tau' `tolerance_coef'
restore

// ============================================================================
// Case 3: Sharp-Frechet, nq=20
// ============================================================================
di _n as text "{hline 70}"
di as text "Case 3: Sharp-Frechet, nq=20, Epanechnikov"
di as text "{hline 70}"

import delimited "ref_data_sharp_simple_20.csv", clear
unab qvars : Q*
local nq : word count `qvars'

local quantiles
forvalues i = 1/`nq' {
    local q = `i' / (`nq' + 1)
    local quantiles `quantiles' `q'
}

r3d X `qvars', cutoff(0) method(frechet) polynomial(2) kernel(epanechnikov) ///
    quantiles(`quantiles') nograph

preserve
import delimited "ref_results_sharp_frechet_20.csv", clear
local max_diff_tau = 0
forvalues i = 1/`nq' {
    local ref_tau = tau[`i']
    local stata_tau = _b[q`i']
    local diff = abs(`stata_tau' - `ref_tau')
    if `diff' > `max_diff_tau' local max_diff_tau = `diff'
}
local pass = (`max_diff_tau' < `tolerance_coef')
report_test "Sharp-Frechet tau" `pass' `max_diff_tau' `tolerance_coef'
restore

// ============================================================================
// Case 5: Fuzzy-Simple, nq=20
// ============================================================================
di _n as text "{hline 70}"
di as text "Case 5: Fuzzy-Simple, nq=20, Epanechnikov"
di as text "{hline 70}"

import delimited "ref_data_fuzzy_simple_20.csv", clear
unab qvars : Q*
local nq : word count `qvars'

local quantiles
forvalues i = 1/`nq' {
    local q = `i' / (`nq' + 1)
    local quantiles `quantiles' `q'
}

r3d X `qvars', cutoff(0) method(simple) polynomial(2) kernel(epanechnikov) ///
    quantiles(`quantiles') fuzzy(t_treat) nograph

preserve
import delimited "ref_results_fuzzy_simple_20.csv", clear
local max_diff_tau = 0
forvalues i = 1/`nq' {
    local ref_tau = tau[`i']
    local stata_tau = _b[q`i']
    local diff = abs(`stata_tau' - `ref_tau')
    if `diff' > `max_diff_tau' local max_diff_tau = `diff'
}
local pass = (`max_diff_tau' < `tolerance_coef')
report_test "Fuzzy-Simple tau" `pass' `max_diff_tau' `tolerance_coef'
restore

// ============================================================================
// Summary
// ============================================================================
di _n as text "{hline 70}"
di as text "EQUIVALENCE TEST SUMMARY"
di as text "{hline 70}"
di as text "Passed: " as result ${n_pass} as text " / " as result ${n_total}
if ${n_fail} > 0 {
    di as error "Failed: " as result ${n_fail}
}
else {
    di as text "All tests passed."
}
di as text "{hline 70}"
