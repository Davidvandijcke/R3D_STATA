/*
 * r3d_mata.mata - Core Mata functions for R3D package
 * Implements regression discontinuity with distributional outcomes
 */

version 18.0
mata:
mata set matastrict off

// ============================================================================
// PLUGIN INTERFACE
// ============================================================================

// Forward declarations for plugin interface (defined in r3d_plugin_interface.mata)
// These are external functions that will be compiled separately

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

// Kernel function implementations
real scalar r3d_kernel(real scalar u, real scalar kernel_type)
{
    if (kernel_type == 1) {  // Triangular
        return((abs(u) <= 1) ? 1 - abs(u) : 0)
    }
    else if (kernel_type == 2) {  // Epanechnikov
        return((abs(u) <= 1) ? 0.75 * (1 - u^2) : 0)
    }
    else {  // Uniform
        return((abs(u) <= 1) ? 0.5 : 0)
    }
}

// Compute variance of a vector
real scalar r3d_variance(real colvector x)
{
    real scalar n, mean_x
    n = rows(x)
    if (n <= 1) return(0)
    mean_x = mean(x)
    return(sum((x :- mean_x):^2) / (n - 1))
}

// ============================================================================
// QUANTILE COMPUTATION
// ============================================================================

// Compute empirical quantiles from distribution data
void r3d_compute_quantiles(string scalar yvars, string scalar quantiles_str,
                          string scalar touse, string scalar matname,
                          string scalar weightvar)
{
    real matrix Y, Q
    real colvector q_grid, w
    real scalar i, j, n, nq, ny
    string rowvector varlist

    // Parse variable list
    varlist = tokens(yvars)
    ny = cols(varlist)

    // Get quantile grid
    q_grid = strtoreal(tokens(quantiles_str))'
    nq = rows(q_grid)

    // Get data
    st_view(Y, ., varlist, touse)
    n = rows(Y)

    // Get weights if specified
    if (weightvar != "") {
        st_view(w, ., weightvar, touse)
    }
    else {
        w = J(n, 1, 1)
    }

    // Initialize quantile matrix
    Q = J(n, nq, .)

    // Compute quantiles for each observation
    for (i = 1; i <= n; i++) {
        // Get non-missing values for this observation
        real rowvector yi, wi
        real scalar k

        yi = select(Y[i,], Y[i,] :< .)
        k = cols(yi)

        if (k > 0) {
            // Create weights for this observation's values
            if (weightvar != "") {
                wi = J(1, k, w[i]/k)  // Equal weight within observation
            }
            else {
                wi = J(1, k, 1/k)
            }

            // Sort values and weights together
            real matrix sorted
            sorted = sort((yi' , wi'), 1)

            // Compute weighted quantiles
            real colvector cumw
            cumw = runningsum(sorted[,2])
            cumw = cumw / cumw[k]  // Normalize

            for (j = 1; j <= nq; j++) {
                real scalar idx
                // Find index where cumulative weight exceeds quantile
                for (idx = 1; idx <= k; idx++) {
                    if (cumw[idx] >= q_grid[j]) break
                }
                if (idx > k) idx = k
                Q[i,j] = sorted[idx,1]
            }
        }
    }

    // Store result
    st_matrix(matname, Q)
}

// ============================================================================
// BANDWIDTH SELECTION
// ============================================================================

// Estimate density at cutoff using Silverman's rule with selected kernel
real scalar r3d_estimate_density(real colvector X, real scalar kernel_type)
{
    real scalar n, h, f_hat
    real colvector X_centered

    n = rows(X)
    X_centered = X :- 0  // Already centered at cutoff

    // Silverman's rule for bandwidth
    real scalar iqr, sd, q25, q75
    real colvector X_sorted
    X_sorted = sort(X_centered, 1)
    q25 = X_sorted[ceil(0.25*n)]
    q75 = X_sorted[ceil(0.75*n)]
    iqr = q75 - q25
    sd = sqrt(r3d_variance(X_centered))
    real scalar scale
    scale = sd
    if (iqr > 0) {
        real scalar alt
        alt = iqr/1.349
        if (alt < scale) scale = alt
    }
    if (scale <= 0) scale = sd
    if (scale <= 0) scale = 1
    h = 1.06 * scale * n^(-0.2)

    // Estimate density at zero
    real colvector K
    K = J(n, 1, 0)
    real scalar i
    for (i = 1; i <= n; i++) {
        K[i] = r3d_kernel(X_centered[i] / h, kernel_type)
    }
    f_hat = mean(K) / h

    return(f_hat)
}

// Compute kernel matrices required for bandwidth selection (Appendix A.3)
void r3d_kernel_matrices(real scalar s, real scalar kernel_type,
                         real matrix Gamma_plus, real matrix Gamma_minus,
                         real colvector Lambda_plus, real colvector Lambda_minus,
                         real matrix Psi_plus, real matrix Psi_minus)
{
    real scalar m, step, i, j, dim
    real rowvector grid

    dim = s + 1
    m = 1000
    grid = range(-1, 1, m)
    step = 2.0 / (m - 1)

    Gamma_plus  = J(dim, dim, 0)
    Gamma_minus = J(dim, dim, 0)
    Psi_plus    = J(dim, dim, 0)
    Psi_minus   = J(dim, dim, 0)
    Lambda_plus  = J(dim, 1, 0)
    Lambda_minus = J(dim, 1, 0)

    for (i = 1; i <= m; i++) {
        real scalar u, k_val
        real colvector basis

        u = grid[i]
        k_val = r3d_kernel(u, kernel_type)
        if (k_val <= 0) continue

        basis = J(dim, 1, 1)
        for (j = 2; j <= dim; j++) {
            basis[j] = basis[j-1] * u
        }

        if (u >= 0) {
            Gamma_plus = Gamma_plus + k_val * (basis * basis') * step
            Psi_plus   = Psi_plus   + (k_val^2) * (basis * basis') * step
            Lambda_plus = Lambda_plus + (u^(s+1)) * k_val * basis * step
        }
        else {
            Gamma_minus = Gamma_minus + k_val * (basis * basis') * step
            Psi_minus   = Psi_minus   + (k_val^2) * (basis * basis') * step
            Lambda_minus = Lambda_minus + (u^(s+1)) * k_val * basis * step
        }
    }
}

// Fit global polynomial and extract derivatives
real rowvector r3d_fit_global_poly(real colvector Y, real colvector X, real scalar order,
                                  real scalar resid_var)
{
    real matrix Xpoly
    real colvector coefs
    real scalar i, n

    n = rows(X)

    // Build polynomial design matrix
    Xpoly = J(n, order+1, 1)
    for (i = 2; i <= order+1; i++) {
        Xpoly[,i] = X:^(i-1)
    }

    // Fit via OLS
    real rowvector derivs
    if (hasmissing(Y) | hasmissing(X)) {
        derivs = J(1, order, 0)
        resid_var = 1
        return(derivs)
    }

    real matrix XtX
    XtX = cross(Xpoly, Xpoly)

    if (det(XtX) != 0) {
        coefs = lusolve(XtX, cross(Xpoly, Y))

        // Extract derivatives adjusted by factorial
        derivs = J(1, order, 0)
        for (i = 1; i <= order; i++) {
            if (i+1 <= rows(coefs)) {
                derivs[i] = coefs[i+1] * factorial(i)
            }
        }

        // Compute residual variance
        real colvector resid
        resid = Y - Xpoly * coefs
        resid_var = r3d_variance(resid)
    }
    else {
        derivs = J(1, order, 0)
        resid_var = 1
    }

    return(derivs)
}

// Core bandwidth selection routine mirroring the R implementation
void r3d_bandwidth_select(string scalar xvar, string scalar Qmat,
                          string scalar quantiles_str, string scalar tvar,
                          real scalar p, real scalar s, real scalar kernel_type,
                          string scalar method, string scalar touse,
                          real scalar is_fuzzy, real scalar coverage,
                          string scalar weightvar)
{
    real matrix X, Q, T
    real rowvector q_grid

    st_view(X, ., xvar, touse)
    Q = st_matrix(Qmat)

    if (quantiles_str == "") {
        errprintf("Quantile grid not supplied to bandwidth selector\n")
        error(498)
    }

    q_grid = strtoreal(tokens(quantiles_str))
    if (cols(q_grid) != cols(Q)) {
        errprintf("Quantile grid length does not match Q matrix\n")
        error(498)
    }

    real scalar n, nq
    n = rows(X)
    nq = cols(Q)

    if (is_fuzzy) {
        st_view(T, ., tvar, touse)
    }

    real scalar sigma_X, f_X_hat
    sigma_X = sqrt(r3d_variance(X))
    if (sigma_X <= 0) sigma_X = 1
    f_X_hat = r3d_estimate_density(X, kernel_type)
    if (f_X_hat <= 0) f_X_hat = 1

    real matrix Gamma_plus, Gamma_minus, Psi_plus, Psi_minus
    real colvector Lambda_plus, Lambda_minus
    r3d_kernel_matrices(s, kernel_type, Gamma_plus, Gamma_minus,
        Lambda_plus, Lambda_minus, Psi_plus, Psi_minus)

    real matrix I_s
    I_s = I(s+1)
    if (abs(det(Gamma_plus)) < 1e-10)  Gamma_plus  = Gamma_plus  + 1e-8 * I_s
    if (abs(det(Gamma_minus)) < 1e-10) Gamma_minus = Gamma_minus + 1e-8 * I_s

    real matrix Gamma_plus_inv, Gamma_minus_inv
    Gamma_plus_inv = invsym(Gamma_plus)
    Gamma_minus_inv = invsym(Gamma_minus)

    real colvector idx_plus, idx_minus
    idx_plus = (X :>= 0)
    idx_minus = (X :< 0)

    real colvector pilot_derivs_plus, pilot_derivs_minus
    real colvector pilot_vars_plus, pilot_vars_minus
    pilot_derivs_plus = J(nq, 1, 0)
    pilot_derivs_minus = J(nq, 1, 0)
    pilot_vars_plus = J(nq, 1, 1)
    pilot_vars_minus = J(nq, 1, 1)

    real scalar order_pilot
    order_pilot = s + 1

    real scalar j
    for (j = 1; j <= nq; j++) {
        real colvector Y_plus, X_plus
        Y_plus = select(Q[,j], idx_plus)
        X_plus = select(X, idx_plus)
        if (rows(Y_plus) > order_pilot) {
            real scalar var_plus
            real rowvector derivs_plus
            derivs_plus = r3d_fit_global_poly(Y_plus, X_plus, order_pilot, var_plus)
            if (cols(derivs_plus) >= order_pilot) pilot_derivs_plus[j] = derivs_plus[order_pilot]
            pilot_vars_plus[j] = var_plus
        }

        real colvector Y_minus, X_minus
        Y_minus = select(Q[,j], idx_minus)
        X_minus = select(X, idx_minus)
        if (rows(Y_minus) > order_pilot) {
            real scalar var_minus
            real rowvector derivs_minus
            derivs_minus = r3d_fit_global_poly(Y_minus, X_minus, order_pilot, var_minus)
            if (cols(derivs_minus) >= order_pilot) pilot_derivs_minus[j] = derivs_minus[order_pilot]
            pilot_vars_minus[j] = var_minus
        }
    }

    real scalar pilot_deriv_T_plus, pilot_deriv_T_minus
    real scalar pilot_var_T_plus, pilot_var_T_minus
    pilot_deriv_T_plus = 0
    pilot_deriv_T_minus = 0
    pilot_var_T_plus = 1
    pilot_var_T_minus = 1

    if (is_fuzzy) {
        real colvector T_plus, T_minus
        T_plus = select(T, idx_plus)
        if (rows(T_plus) > order_pilot) {
            real scalar var_T_plus
            real rowvector deriv_T_plus
            deriv_T_plus = r3d_fit_global_poly(T_plus, X_plus, order_pilot, var_T_plus)
            if (cols(deriv_T_plus) >= order_pilot) pilot_deriv_T_plus = deriv_T_plus[order_pilot]
            pilot_var_T_plus = var_T_plus
        }

        T_minus = select(T, idx_minus)
        if (rows(T_minus) > order_pilot) {
            real scalar var_T_minus
            real rowvector deriv_T_minus
            deriv_T_minus = r3d_fit_global_poly(T_minus, X_minus, order_pilot, var_T_minus)
            if (cols(deriv_T_minus) >= order_pilot) pilot_deriv_T_minus = deriv_T_minus[order_pilot]
            pilot_var_T_minus = var_T_minus
        }
    }

    real colvector e0
    e0 = J(s+1, 1, 0)
    e0[1] = 1

    real colvector pilot_h_num
    pilot_h_num = J(nq, 1, 0)

    real scalar factorial_term
    factorial_term = factorial(s + 1)

    for (j = 1; j <= nq; j++) {
        real scalar bias_plus, bias_minus, C_1_0
        real scalar C_1_0_prime

        bias_plus = (e0' * (Gamma_plus_inv * Lambda_plus)) * (pilot_derivs_plus[j] / factorial_term)
        bias_minus = (e0' * (Gamma_minus_inv * Lambda_minus)) * (pilot_derivs_minus[j] / factorial_term)
        C_1_0 = bias_plus - bias_minus

        real colvector temp_plus, temp_minus
        real scalar var_plus_j, var_minus_j
        temp_plus = Gamma_plus_inv * e0
        temp_minus = Gamma_minus_inv * e0

        var_plus_j = pilot_vars_plus[j] * (temp_plus' * (Psi_plus * temp_plus))
        var_minus_j = pilot_vars_minus[j] * (temp_minus' * (Psi_minus * temp_minus))
        C_1_0_prime = (var_plus_j + var_minus_j) / f_X_hat

        if (abs(C_1_0) > 1e-14 && C_1_0_prime > 0) {
            real scalar ratio
            ratio = C_1_0_prime / (2 * (s + 1) * (C_1_0^2))
            pilot_h_num[j] = ratio^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
        }

        if (pilot_h_num[j] <= 0) {
            pilot_h_num[j] = 1.06 * sigma_X * n^(-1 / (2 * s + 1))
        }
    }

    real scalar pilot_h_den
    pilot_h_den = .
    if (is_fuzzy) {
        real scalar bias_T_plus, bias_T_minus, C_T_0
        real scalar C_T_0_prime

        bias_T_plus = (e0' * (Gamma_plus_inv * Lambda_plus)) * (pilot_deriv_T_plus / factorial_term)
        bias_T_minus = (e0' * (Gamma_minus_inv * Lambda_minus)) * (pilot_deriv_T_minus / factorial_term)
        C_T_0 = bias_T_plus - bias_T_minus

        real colvector temp_plus_T, temp_minus_T
        temp_plus_T = Gamma_plus_inv * e0
        temp_minus_T = Gamma_minus_inv * e0

        real scalar var_T_plus_calc, var_T_minus_calc
        var_T_plus_calc = pilot_var_T_plus * (temp_plus_T' * (Psi_plus * temp_plus_T))
        var_T_minus_calc = pilot_var_T_minus * (temp_minus_T' * (Psi_minus * temp_minus_T))
        C_T_0_prime = (var_T_plus_calc + var_T_minus_calc) / f_X_hat

        if (abs(C_T_0) > 1e-14 && C_T_0_prime > 0) {
            real scalar ratio_T
            ratio_T = C_T_0_prime / (2 * (s + 1) * (C_T_0^2))
            pilot_h_den = ratio_T^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
        }

        if (!(pilot_h_den > 0)) {
            pilot_h_den = 1.06 * sigma_X * n^(-1 / (2 * s + 1))
        }
    }

    real matrix alpha_plus_pilot, alpha_minus_pilot, w_plus_pilot, w_minus_pilot
    real matrix h_pilot_vec
    h_pilot_vec = pilot_h_num

    if (method == "frechet") {
        real scalar h_scalar
        h_scalar = mean(h_pilot_vec)
        if (!(h_scalar > 0)) h_scalar = 1.06 * sigma_X * n^(-1 / (2 * s + 1))
        h_pilot_vec = J(nq, 1, h_scalar)
    }

    if (r3d_plugin_available()) {
        real scalar rc_plus, rc_minus
        rc_plus = r3d_plugin_locweights(X, Q, h_pilot_vec, s, kernel_type, 1, alpha_plus_pilot, w_plus_pilot)
        rc_minus = r3d_plugin_locweights(X, Q, h_pilot_vec, s, kernel_type, 0, alpha_minus_pilot, w_minus_pilot)
        if (rc_plus != 0 | rc_minus != 0) {
            st_global("R3D_USE_PLUGIN", "0")
            st_global("R3D_PLUGIN_NAME", "")
            r3d_locpoly(X, Q, h_pilot_vec, s, kernel_type, 1, alpha_plus_pilot, w_plus_pilot)
            r3d_locpoly(X, Q, h_pilot_vec, s, kernel_type, 0, alpha_minus_pilot, w_minus_pilot)
        }
    }
    else {
        r3d_locpoly(X, Q, h_pilot_vec, s, kernel_type, 1, alpha_plus_pilot, w_plus_pilot)
        r3d_locpoly(X, Q, h_pilot_vec, s, kernel_type, 0, alpha_minus_pilot, w_minus_pilot)
    }

    real matrix alphaT_plus_pilot, alphaT_minus_pilot, wT_plus_pilot, wT_minus_pilot
    if (is_fuzzy) {
        real matrix h_den_vec
        h_den_vec = J(1, 1, pilot_h_den)

        if (r3d_plugin_available()) {
            real scalar rc_tplus, rc_tminus
            rc_tplus = r3d_plugin_locweights(X, T, h_den_vec, s, kernel_type, 1, alphaT_plus_pilot, wT_plus_pilot)
            rc_tminus = r3d_plugin_locweights(X, T, h_den_vec, s, kernel_type, 0, alphaT_minus_pilot, wT_minus_pilot)
            if (rc_tplus != 0 | rc_tminus != 0) {
                st_global("R3D_USE_PLUGIN", "0")
                st_global("R3D_PLUGIN_NAME", "")
                r3d_locpoly(X, T, h_den_vec, s, kernel_type, 1, alphaT_plus_pilot, wT_plus_pilot)
                r3d_locpoly(X, T, h_den_vec, s, kernel_type, 0, alphaT_minus_pilot, wT_minus_pilot)
            }
        }
        else {
            r3d_locpoly(X, T, h_den_vec, s, kernel_type, 1, alphaT_plus_pilot, wT_plus_pilot)
            r3d_locpoly(X, T, h_den_vec, s, kernel_type, 0, alphaT_minus_pilot, wT_minus_pilot)
        }
    }

    real colvector B_plus, B_minus, V_plus, V_minus
    B_plus = J(nq, 1, 0)
    B_minus = J(nq, 1, 0)
    V_plus = J(nq, 1, 0)
    V_minus = J(nq, 1, 0)

    for (j = 1; j <= nq; j++) {
        B_plus[j] = alpha_plus_pilot[s+1, j] / factorial_term
        B_minus[j] = alpha_minus_pilot[s+1, j] / factorial_term

        real scalar h_use, weight_sum, k
        h_use = h_pilot_vec[j,1]
        if (!(h_use > 0)) continue

        real colvector idxp, idxm
        idxp = (w_plus_pilot[,j] :> 0)
        idxm = (w_minus_pilot[,j] :> 0)

        if (sum(idxp) > 0) {
            real colvector Xscaled_p
            real matrix basis_p
            real colvector fitted_p, resid_p

            Xscaled_p = select(X, idxp) / h_use
            basis_p = J(rows(Xscaled_p), s+1, 1)
            for (k = 2; k <= s+1; k++) basis_p[,k] = basis_p[,k-1] :* Xscaled_p
            fitted_p = basis_p * alpha_plus_pilot[,j]
            resid_p = select(Q[,j], idxp) - fitted_p
            weight_sum = sum(select(w_plus_pilot[,j], idxp))
            if (weight_sum > 0) V_plus[j] = sum(select(w_plus_pilot[,j], idxp) :* (resid_p:^2)) / weight_sum
        }

        if (sum(idxm) > 0) {
            real colvector Xscaled_m
            real matrix basis_m
            real colvector fitted_m, resid_m
            Xscaled_m = select(X, idxm) / h_use
            basis_m = J(rows(Xscaled_m), s+1, 1)
            for (k = 2; k <= s+1; k++) basis_m[,k] = basis_m[,k-1] :* Xscaled_m
            fitted_m = basis_m * alpha_minus_pilot[,j]
            resid_m = select(Q[,j], idxm) - fitted_m
            weight_sum = sum(select(w_minus_pilot[,j], idxm))
            if (weight_sum > 0) V_minus[j] = sum(select(w_minus_pilot[,j], idxm) :* (resid_m:^2)) / weight_sum
        }
    }

    real scalar B_plus_den, B_minus_den, V_plus_den, V_minus_den
    B_plus_den = 0
    B_minus_den = 0
    V_plus_den = 0
    V_minus_den = 0

    if (is_fuzzy) {
        B_plus_den = alphaT_plus_pilot[s+1,1] / factorial_term
        B_minus_den = alphaT_minus_pilot[s+1,1] / factorial_term

        real colvector idxp_den, idxm_den
        real scalar weight_sum_den, k_den
        idxp_den = (wT_plus_pilot :> 0)
        idxm_den = (wT_minus_pilot :> 0)

        if (sum(idxp_den) > 0) {
            real colvector Xscaled_p
            real matrix basis_p
            real colvector fitted_p, resid_p

            Xscaled_p = select(X, idxp_den) / pilot_h_den
            basis_p = J(rows(Xscaled_p), s+1, 1)
            for (k_den = 2; k_den <= s+1; k_den++) basis_p[,k_den] = basis_p[,k_den-1] :* Xscaled_p
            fitted_p = basis_p * alphaT_plus_pilot
            resid_p = select(T, idxp_den) - fitted_p
            weight_sum_den = sum(select(wT_plus_pilot, idxp_den))
            if (weight_sum_den > 0) V_plus_den = sum(select(wT_plus_pilot, idxp_den) :* (resid_p:^2)) / weight_sum_den
        }

        if (sum(idxm_den) > 0) {
            real colvector Xscaled_m
            real matrix basis_m
            real colvector fitted_m, resid_m

            Xscaled_m = select(X, idxm_den) / pilot_h_den
            basis_m = J(rows(Xscaled_m), s+1, 1)
            for (k_den = 2; k_den <= s+1; k_den++) basis_m[,k_den] = basis_m[,k_den-1] :* Xscaled_m
            fitted_m = basis_m * alphaT_minus_pilot
            resid_m = select(T, idxm_den) - fitted_m
            weight_sum_den = sum(select(wT_minus_pilot, idxm_den))
            if (weight_sum_den > 0) V_minus_den = sum(select(wT_minus_pilot, idxm_den) :* (resid_m:^2)) / weight_sum_den
        }
    }

    real colvector h_star_num
    h_star_num = J(nq, 1, 0)

    if (method == "simple") {
        for (j = 1; j <= nq; j++) {
            real scalar bias_diff, var_sum
            bias_diff = B_plus[j] - B_minus[j]
            var_sum = (V_plus[j] + V_minus[j]) / f_X_hat
            if (abs(bias_diff) > 1e-14 && var_sum > 0) {
                real scalar ratio
                ratio = var_sum / (2 * (s + 1) * (bias_diff^2))
                h_star_num[j] = ratio^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
            }
            if (h_star_num[j] <= 0) h_star_num[j] = 1.06 * sigma_X * n^(-1 / (2 * s + 1))
        }
    }
    else {
        real scalar dq, A_s, B_s
        dq = (max(q_grid) - min(q_grid)) / (nq - 1)
        if (dq <= 0) dq = 1.0 / nq
        A_s = sum((B_plus - B_minus):^2) * dq
        B_s = sum((V_plus + V_minus)) * dq / f_X_hat
        real scalar h_val
        if (A_s > 1e-14 && B_s > 0) {
            real scalar ratio
            ratio = B_s / (2 * (s + 1) * A_s)
            h_val = ratio^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
        }
        else {
            h_val = 1.06 * sigma_X * n^(-1 / (2 * s + 1))
        }
        h_star_num = J(nq, 1, h_val)
    }

    real scalar h_star_den
    h_star_den = .
    if (is_fuzzy) {
        real scalar bias_diff_den, var_sum_den
        bias_diff_den = B_plus_den - B_minus_den
        var_sum_den = (V_plus_den + V_minus_den) / f_X_hat
        if (abs(bias_diff_den) > 1e-14 && var_sum_den > 0) {
            real scalar ratio_den
            ratio_den = var_sum_den / (2 * (s + 1) * (bias_diff_den^2))
            h_star_den = ratio_den^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
        }
        if (!(h_star_den > 0)) h_star_den = 1.06 * sigma_X * n^(-1 / (2 * s + 1))
    }

    if (coverage) {
        real scalar shrink
        shrink = n^(-s / ((2 * s + 3) * (s + 3)))
        h_star_num = h_star_num :* shrink
        if (is_fuzzy) h_star_den = h_star_den * shrink
    }

    st_matrix("r(h_num)", h_star_num)
    st_matrix("r(pilot_num)", pilot_h_num)
    st_matrix("r(B_plus)", B_plus)
    st_matrix("r(B_minus)", B_minus)
    st_matrix("r(V_plus)", V_plus)
    st_matrix("r(V_minus)", V_minus)
    st_numscalar("r(f_X_hat)", f_X_hat)

    if (is_fuzzy) {
        st_numscalar("r(h_den)", h_star_den)
        st_numscalar("r(pilot_den)", pilot_h_den)
        st_numscalar("r(B_plus_den)", B_plus_den)
        st_numscalar("r(B_minus_den)", B_minus_den)
        st_numscalar("r(V_plus_den)", V_plus_den)
        st_numscalar("r(V_minus_den)", V_minus_den)
    }
}

// Validate and broadcast user-supplied bandwidth vector
real scalar r3d_prepare_bandwidth_matrix(string scalar matname, real scalar nq,
                                        string scalar method)
{
    real matrix H
    H = st_matrix(matname)
    if (rows(H) == 1 & cols(H) == 1) {
        H = J(nq, 1, H[1,1])
    }
    else if (rows(H) == 1 & cols(H) == nq) {
        H = H'
    }
    else if (!(rows(H) == nq & cols(H) == 1)) {
        return(198)
    }

    if (method == "frechet") {
        real scalar unique
        unique = rows( uniqrows(H) )
        if (unique > 1) {
            return(198)
        }
        H = J(nq, 1, H[1,1])
    }

    st_matrix(matname, H)
    return(0)
}

// ============================================================================
// LOCAL POLYNOMIAL REGRESSION
// ============================================================================

// Core local polynomial regression
void r3d_locpoly(real colvector X, real matrix Y, real matrix h_mat,
                real scalar p, real scalar kernel_type, real scalar side,
                real matrix alpha, real matrix weights)
{
    real scalar n, nq, i, j, k
    real matrix XWX, XWY
    real colvector basis
    real scalar u, w_i

    n = rows(X)
    nq = cols(Y)

    // Normalise bandwidth input to column vector of length nq
    if (rows(h_mat) == 1 & cols(h_mat) == 1) {
        h_mat = J(nq, 1, h_mat[1,1])
    }
    else if (rows(h_mat) == 1 & cols(h_mat) == nq) {
        h_mat = h_mat'
    }
    else if (!(rows(h_mat) == nq & cols(h_mat) == 1)) {
        stata("display as error \"r3d_locpoly(): bandwidth vector has incorrect dimensions\"")
        return
    }

    // Initialize outputs
    alpha = J(p+1, nq, 0)
    weights = J(n, nq, 0)

    for (j = 1; j <= nq; j++) {
        real scalar h_j
        h_j = h_mat[j,1]
        if (h_j <= 0) continue

        // Initialize matrices
        XWX = J(p+1, p+1, 0)
        XWY = J(p+1, 1, 0)

        // Compute weighted matrices
        for (i = 1; i <= n; i++) {
            // Skip observations on wrong side
            if (side == 1 & X[i] < 0) continue
            if (side == 0 & X[i] >= 0) continue

            // Skip if Y is missing
            if (Y[i,j] >= .) continue

            // Compute kernel weight
            u = X[i] / h_j
            w_i = r3d_kernel(u, kernel_type)

            if (w_i > 0) {
                // Store weight
                weights[i,j] = w_i

                // Compute polynomial basis
                basis = J(p+1, 1, 1)
                for (k = 2; k <= p+1; k++) {
                    basis[k] = basis[k-1] * (X[i]/h_j)
                }

                // Accumulate
                XWX = XWX + w_i * basis * basis'
                XWY = XWY + w_i * basis * Y[i,j]
            }
        }

        // Solve for coefficients
        if (det(XWX) != 0) {
            alpha[,j] = lusolve(XWX, XWY)
        }
        else {
            alpha[,j] = J(p+1, 1, 0)
        }
    }
}

// ============================================================================
// ISOTONIC REGRESSION
// ============================================================================

// Pool adjacent violators algorithm for isotonic regression
real colvector r3d_isotonic(real colvector y)
{
    real scalar n, i
    real colvector y_iso

    n = rows(y)
    if (n <= 1) return(y)

    y_iso = y

    // Pool adjacent violators algorithm
    i = 1
    while (i < n) {
        if (y_iso[i] > y_iso[i+1]) {
            // Find violating block
            real scalar j, sum_y, len
            j = i + 1
            sum_y = y_iso[i] + y_iso[i+1]
            len = 2

            // Extend block while violations continue
            while (j < n & sum_y/len > y_iso[j+1]) {
                j++
                sum_y = sum_y + y_iso[j]
                len++
            }

            // Pool the block
            real scalar avg, k
            avg = sum_y / len
            for (k = i; k <= j; k++) {
                y_iso[k] = avg
            }

            // Backtrack to check previous values
            if (i > 1) {
                i--
            }
        }
        else {
            i++
        }
    }

    return(y_iso)
}

// ============================================================================
// MAIN ESTIMATION FUNCTIONS
// ============================================================================

// Core estimator covering both simple and Frechet methods
void r3d_estimate_core(string scalar method,
                      string scalar xvar, string scalar Qmat, string scalar tvar,
                      string scalar hnum_matname, real scalar h_den,
                      real scalar p, real scalar kernel_type,
                      string scalar tau_name, string scalar se_name,
                      string scalar alpha_plus_name, string scalar alpha_minus_name,
                      string scalar w_plus_name, string scalar w_minus_name,
                      string scalar e1_name, string scalar e2_name,
                      string scalar int_plus_name, string scalar int_minus_name,
                      string scalar alpha_t_plus_name, string scalar alpha_t_minus_name,
                      string scalar w_t_plus_name, string scalar w_t_minus_name,
                      string scalar denom_name,
                      string scalar touse, real scalar is_fuzzy)
{
    real matrix X, Q
    st_view(X, ., xvar, touse)
    Q = st_matrix(Qmat)

    real scalar n, nq
    n = rows(X)
    nq = cols(Q)

    real matrix h_num_vec
    h_num_vec = st_matrix(hnum_matname)

    if (rows(h_num_vec) == 1 & cols(h_num_vec) == 1) {
        h_num_vec = J(nq, 1, h_num_vec[1,1])
    }
    else if (rows(h_num_vec) == 1 & cols(h_num_vec) == nq) {
        h_num_vec = h_num_vec'
    }

    if (rows(h_num_vec) != nq | cols(h_num_vec) != 1) {
        stata("display as error \"Bandwidth vector dimension mismatch in r3d_estimate_core()\"")
        return
    }

    real matrix alpha_plus, alpha_minus, w_plus, w_minus

    if (r3d_plugin_available()) {
        real scalar rc_plus, rc_minus
        rc_plus = r3d_plugin_locweights(X, Q, h_num_vec, p, kernel_type, 1, alpha_plus, w_plus)
        rc_minus = r3d_plugin_locweights(X, Q, h_num_vec, p, kernel_type, 0, alpha_minus, w_minus)
        if (rc_plus != 0 | rc_minus != 0) {
            st_global("R3D_USE_PLUGIN", "0")
            st_global("R3D_PLUGIN_NAME", "")
            r3d_locpoly(X, Q, h_num_vec, p, kernel_type, 1, alpha_plus, w_plus)
            r3d_locpoly(X, Q, h_num_vec, p, kernel_type, 0, alpha_minus, w_minus)
        }
    }
    else {
        r3d_locpoly(X, Q, h_num_vec, p, kernel_type, 1, alpha_plus, w_plus)
        r3d_locpoly(X, Q, h_num_vec, p, kernel_type, 0, alpha_minus, w_minus)
    }

    real matrix T
    real matrix alpha_plus_t, alpha_minus_t, w_plus_t, w_minus_t
    real matrix e2
    real scalar denominator
    denominator = 1
    e2 = J(n, nq, 0)

    if (is_fuzzy) {
        st_view(T, ., tvar, touse)

        real matrix h_den_vec
        h_den_vec = J(1, 1, h_den)

        if (r3d_plugin_available()) {
            real scalar rc_tplus, rc_tminus
            rc_tplus = r3d_plugin_locweights(X, T, h_den_vec, p, kernel_type, 1, alpha_plus_t, w_plus_t)
            rc_tminus = r3d_plugin_locweights(X, T, h_den_vec, p, kernel_type, 0, alpha_minus_t, w_minus_t)
            if (rc_tplus != 0 | rc_tminus != 0) {
                st_global("R3D_USE_PLUGIN", "0")
                st_global("R3D_PLUGIN_NAME", "")
                r3d_locpoly(X, T, h_den_vec, p, kernel_type, 1, alpha_plus_t, w_plus_t)
                r3d_locpoly(X, T, h_den_vec, p, kernel_type, 0, alpha_minus_t, w_minus_t)
            }
        }
        else {
            r3d_locpoly(X, T, h_den_vec, p, kernel_type, 1, alpha_plus_t, w_plus_t)
            r3d_locpoly(X, T, h_den_vec, p, kernel_type, 0, alpha_minus_t, w_minus_t)
        }

        denominator = alpha_plus_t[1,1] - alpha_minus_t[1,1]
    }

    // Predicted values on each side
    real matrix Eplus, Eminus, Eplus_final, Eminus_final
    Eplus = J(n, nq, .)
    Eminus = J(n, nq, .)
    real scalar j, k

    for (j = 1; j <= nq; j++) {
        real scalar h_j
        h_j = h_num_vec[j,1]
        if (!(h_j > 0)) continue

        real colvector idxp, idxm
        idxp = (w_plus[,j] :> 0)
        idxm = (w_minus[,j] :> 0)

        if (sum(idxp) > 0) {
            real colvector Xscaled_p
            real matrix basis_p
            real colvector fitted_p
            Xscaled_p = select(X, idxp) / h_j
            basis_p = J(rows(Xscaled_p), p+1, 1)
            for (k = 2; k <= p+1; k++) {
                basis_p[,k] = basis_p[,k-1] :* Xscaled_p
            }
            fitted_p = basis_p * alpha_plus[,j]
            Eplus[idxp, j] = fitted_p
        }

        if (sum(idxm) > 0) {
            real colvector Xscaled_m
            real matrix basis_m
            real colvector fitted_m
            Xscaled_m = select(X, idxm) / h_j
            basis_m = J(rows(Xscaled_m), p+1, 1)
            for (k = 2; k <= p+1; k++) {
                basis_m[,k] = basis_m[,k-1] :* Xscaled_m
            }
            fitted_m = basis_m * alpha_minus[,j]
            Eminus[idxm, j] = fitted_m
        }
    }

    Eplus_final = Eplus
    Eminus_final = Eminus

    if (method == "frechet") {
        real scalar i
        for (i = 1; i <= n; i++) {
            if (sum(w_plus[i,] :> 0)) {
                real colvector rowfit
                rowfit = Eplus[i,]'
                rowfit = r3d_isotonic(rowfit)
                Eplus_final[i,] = rowfit'
            }
            if (sum(w_minus[i,] :> 0)) {
                real colvector rowfit
                rowfit = Eminus[i,]'
                rowfit = r3d_isotonic(rowfit)
                Eminus_final[i,] = rowfit'
            }
        }
    }

    // Residuals for outcome
    real matrix e1
    e1 = J(n, nq, 0)
    for (j = 1; j <= nq; j++) {
        real colvector idxp, idxm
        idxp = (w_plus[,j] :> 0)
        idxm = (w_minus[,j] :> 0)
        if (sum(idxp) > 0) {
            e1[idxp, j] = select(Q[,j], idxp) - Eplus_final[idxp, j]
        }
        if (sum(idxm) > 0) {
            e1[idxm, j] = select(Q[,j], idxm) - Eminus_final[idxm, j]
        }
    }

    // Residuals for treatment (replicated across quantiles)
    if (is_fuzzy) {
        real colvector idxp, idxm
        idxp = (w_plus_t :> 0)
        idxm = (w_minus_t :> 0)
        if (sum(idxp) > 0) {
            real colvector Xscaled_p
            real matrix basis_p
            real colvector fitted_p, resid_p
            Xscaled_p = select(X, idxp) / h_den
            basis_p = J(rows(Xscaled_p), p+1, 1)
            for (k = 2; k <= p+1; k++) {
                basis_p[,k] = basis_p[,k-1] :* Xscaled_p
            }
            fitted_p = basis_p * alpha_plus_t
            resid_p = select(T, idxp) - fitted_p
            e2[idxp,] = resid_p * J(1, nq, 1)
        }
        if (sum(idxm) > 0) {
            real colvector Xscaled_m
            real matrix basis_m
            real colvector fitted_m, resid_m
            Xscaled_m = select(X, idxm) / h_den
            basis_m = J(rows(Xscaled_m), p+1, 1)
            for (k = 2; k <= p+1; k++) {
                basis_m[,k] = basis_m[,k-1] :* Xscaled_m
            }
            fitted_m = basis_m * alpha_minus_t
            resid_m = select(T, idxm) - fitted_m
            e2[idxm,] = resid_m * J(1, nq, 1)
        }
    }

    // Intercepts and rearrangement/isotonic adjustments
    real rowvector int_plus, int_minus
    int_plus = alpha_plus[1,]
    int_minus = alpha_minus[1,]

    if (method == "simple") {
        int_plus = sort(int_plus', 1)'
        int_minus = sort(int_minus', 1)'
    }
    else {
        int_plus = r3d_isotonic(int_plus')'
        int_minus = r3d_isotonic(int_minus')'
    }
    alpha_plus[1,] = int_plus
    alpha_minus[1,] = int_minus

    // Quantile treatment effect
    real rowvector tau
    tau = int_plus - int_minus

    if (is_fuzzy) {
        if (abs(denominator) > 1e-8) {
            tau = tau / denominator
        }
        else {
            tau = J(1, nq, .)
        }
    }

    // Standard errors (approximate)
    real rowvector se
    se = J(1, nq, .)
    for (j = 1; j <= nq; j++) {
        real colvector idxp, idxm
        idxp = (w_plus[,j] :> 0)
        idxm = (w_minus[,j] :> 0)
        real scalar n_plus, n_minus
        n_plus = sum(idxp)
        n_minus = sum(idxm)
        if (n_plus > p+1 & n_minus > p+1) {
            real scalar var_plus, var_minus
            var_plus = r3d_variance(select(e1[,j], idxp))
            var_minus = r3d_variance(select(e1[,j], idxm))
            if (var_plus >= 0 & var_minus >= 0) {
                se[j] = sqrt(var_plus/n_plus + var_minus/n_minus)
                if (is_fuzzy & abs(denominator) > 1e-8) {
                    se[j] = se[j] / abs(denominator)
                }
            }
        }
    }

    // Store results
    st_matrix(tau_name, tau)
    st_matrix(se_name, se)
    st_matrix(alpha_plus_name, alpha_plus)
    st_matrix(alpha_minus_name, alpha_minus)
    st_matrix(w_plus_name, w_plus)
    st_matrix(w_minus_name, w_minus)
    st_matrix(e1_name, e1)
    st_matrix(int_plus_name, int_plus)
    st_matrix(int_minus_name, int_minus)
    st_numscalar(denom_name, denominator)

    if (is_fuzzy) {
        st_matrix(e2_name, e2)
        if (alpha_t_plus_name != "") st_matrix(alpha_t_plus_name, alpha_plus_t)
        if (alpha_t_minus_name != "") st_matrix(alpha_t_minus_name, alpha_minus_t)
        if (w_t_plus_name != "") st_matrix(w_t_plus_name, w_plus_t)
        if (w_t_minus_name != "") st_matrix(w_t_minus_name, w_minus_t)
    }
    else {
        st_matrix(e2_name, J(n, nq, 0))
        if (alpha_t_plus_name != "") st_matrix(alpha_t_plus_name, J(p+1, 1, .))
        if (alpha_t_minus_name != "") st_matrix(alpha_t_minus_name, J(p+1, 1, .))
        if (w_t_plus_name != "") st_matrix(w_t_plus_name, J(n, 1, 0))
        if (w_t_minus_name != "") st_matrix(w_t_minus_name, J(n, 1, 0))
    }
}

// Simple method wrapper
void r3d_simple(string scalar xvar, string scalar Qmat, string scalar tvar,
               string scalar hnum_matname, real scalar h_den, real scalar p,
               real scalar kernel_type, string scalar tau_name, string scalar se_name,
               string scalar alpha_plus_name, string scalar alpha_minus_name,
               string scalar w_plus_name, string scalar w_minus_name,
               string scalar e1_name, string scalar e2_name,
               string scalar int_plus_name, string scalar int_minus_name,
               string scalar alpha_t_plus_name, string scalar alpha_t_minus_name,
               string scalar w_t_plus_name, string scalar w_t_minus_name,
               string scalar denom_name,
               string scalar touse, real scalar is_fuzzy)
{
    r3d_estimate_core("simple", xvar, Qmat, tvar, hnum_matname, h_den, p, kernel_type,
        tau_name, se_name, alpha_plus_name, alpha_minus_name,
        w_plus_name, w_minus_name, e1_name, e2_name,
        int_plus_name, int_minus_name,
        alpha_t_plus_name, alpha_t_minus_name,
        w_t_plus_name, w_t_minus_name,
        denom_name,
        touse, is_fuzzy);
}

// Frechet method wrapper
void r3d_frechet(string scalar xvar, string scalar Qmat, string scalar tvar,
                string scalar hnum_matname, real scalar h_den, real scalar p,
                real scalar kernel_type, string scalar tau_name, string scalar se_name,
                string scalar alpha_plus_name, string scalar alpha_minus_name,
                string scalar w_plus_name, string scalar w_minus_name,
                string scalar e1_name, string scalar e2_name,
                string scalar int_plus_name, string scalar int_minus_name,
                string scalar alpha_t_plus_name, string scalar alpha_t_minus_name,
                string scalar w_t_plus_name, string scalar w_t_minus_name,
                string scalar denom_name,
                string scalar touse, real scalar is_fuzzy)
{
    r3d_estimate_core("frechet", xvar, Qmat, tvar, hnum_matname, h_den, p, kernel_type,
        tau_name, se_name, alpha_plus_name, alpha_minus_name,
        w_plus_name, w_minus_name, e1_name, e2_name,
        int_plus_name, int_minus_name,
        alpha_t_plus_name, alpha_t_minus_name,
        w_t_plus_name, w_t_minus_name,
        denom_name,
        touse, is_fuzzy);
}

// ============================================================================
// BOOTSTRAP INFERENCE
// ============================================================================

// Multiplier bootstrap for uniform confidence bands and hypothesis tests
void r3d_bootstrap(string scalar xvar, string scalar Qmat, string scalar tvar,
                  string scalar hnum_matname, real scalar h_den, real scalar p,
                  real scalar kernel_type, real scalar B, real scalar level,
                  string scalar alpha_plus_name, string scalar alpha_minus_name,
                  string scalar w_plus_name, string scalar w_minus_name,
                  string scalar e1_name, string scalar e2_name,
                  string scalar int_plus_name, string scalar int_minus_name,
                  string scalar tau_name, string scalar denom_name,
                  real scalar method_code, string scalar quantiles_str,
                  string scalar tests_str, string scalar test_ranges_str,
                  string scalar bs_results_name, string scalar pvals_name,
                  string scalar touse, real scalar is_fuzzy)
{
    real matrix X
    st_view(X, ., xvar, touse)
    real scalar n
    n = rows(X)

    real matrix h_num_vec
    h_num_vec = st_matrix(hnum_matname)
    if (rows(h_num_vec) == 1 & cols(h_num_vec) == 1) {
        h_num_vec = J(cols(st_matrix(tau_name)), 1, h_num_vec[1,1])
    }
    else if (rows(h_num_vec) == 1 & cols(h_num_vec) > 1) {
        h_num_vec = h_num_vec'
    }

    real matrix alpha_plus, alpha_minus, w_plus, w_minus
    alpha_plus = st_matrix(alpha_plus_name)
    alpha_minus = st_matrix(alpha_minus_name)
    w_plus = st_matrix(w_plus_name)
    w_minus = st_matrix(w_minus_name)

    real matrix e1_mat
    e1_mat = st_matrix(e1_name)

    real rowvector int_plus, int_minus, tau_orig
    int_plus = st_matrix(int_plus_name)
    int_minus = st_matrix(int_minus_name)
    tau_orig = st_matrix(tau_name)

    real scalar denominator
    denominator = st_numscalar(denom_name)

    real scalar nq
    nq = cols(alpha_plus)

    real matrix e1_w_plus, e1_w_minus
    e1_w_plus = e1_mat :* w_plus
    e1_w_minus = e1_mat :* w_minus

    real matrix e2_w_plus, e2_w_minus
    if (is_fuzzy) {
        real matrix e2_mat
        e2_mat = st_matrix(e2_name)
        e2_w_plus = e2_mat :* w_plus
        e2_w_minus = e2_mat :* w_minus
    }

    real rowvector num_diff
    num_diff = int_plus - int_minus

    real scalar is_vector_h_num
    is_vector_h_num = (method_code == 0)

    // Parse tests
    tests_str = strlower(strtrim(subinstr(tests_str, ",", " ", .)))
    string rowvector tests_tokens
    tests_tokens = tokens(tests_str)

    string rowvector filtered
    filtered = J(1, 0, "")
    for (real scalar j = 1; j <= cols(tests_tokens); j++) {
        if (tests_tokens[j] != "" & tests_tokens[j] != "none") {
            filtered = filtered , tests_tokens[j]
        }
    }
    if (cols(filtered) > 0) {
        tests_tokens = filtered
    }
    else {
        tests_tokens = J(1, 0, "")
    }

    real scalar do_nullity, do_homogeneity, do_gini
    do_nullity = 0
    do_homogeneity = 0
    do_gini = 0
    if (cols(tests_tokens) > 0) {
        do_nullity = any(tests_tokens :== "nullity")
        do_homogeneity = any(tests_tokens :== "homogeneity")
        do_gini = any(tests_tokens :== "gini")
    }

    real matrix tau_boot
    tau_boot = J(B, nq, 0)

    real matrix nu_plus_store, nu_minus_store
    if (do_gini) {
        nu_plus_store = J(B, nq, 0)
        nu_minus_store = J(B, nq, 0)
    }

    real colvector xi
    real rowvector plus_sums, minus_sums, out_sharp

    for (real scalar b = 1; b <= B; b++) {
        xi = rnormal(n, 1, 0, 1)
        plus_sums = xi' * e1_w_plus
        minus_sums = xi' * e1_w_minus
        out_sharp = plus_sums - minus_sums

        if (do_gini) {
            nu_plus_store[b,] = plus_sums
            nu_minus_store[b,] = minus_sums
        }

        if (!is_fuzzy) {
            tau_boot[b,] = out_sharp
        }
        else {
            real scalar denom_term
            denom_term = (xi' * e2_w_plus[,1]) - (xi' * e2_w_minus[,1])
            if (abs(denominator) > 1e-8) {
                real rowvector top
                top = denominator * out_sharp - num_diff :* denom_term
                tau_boot[b,] = top / (denominator^2)
            }
            else {
                tau_boot[b,] = J(1, nq, 0)
            }
        }
    }

    // Uniform confidence bands
    real colvector supvals
    supvals = J(B, 1, 0)
    for (real scalar b = 1; b <= B; b++) {
        supvals[b] = max(abs(tau_boot[b,]))
    }

    real colvector sup_sorted
    sup_sorted = sort(supvals, 1)

    real scalar prob, alpha_level, index_c
    prob = level / 100
    alpha_level = 1 - prob
    index_c = ceil(prob * (B + 1))
    if (index_c < 1) index_c = 1
    if (index_c > B) index_c = B

    real scalar cval
    cval = sup_sorted[index_c]

    real rowvector cb_lower, cb_upper
    cb_lower = tau_orig - cval
    cb_upper = tau_orig + cval

    // Prepare quantile grid and test ranges
    real rowvector q_grid
    q_grid = strtoreal(tokens(quantiles_str))

    real matrix range_pairs
    real rowvector range_tokens
    range_tokens = strtoreal(tokens(test_ranges_str))
    if (cols(range_tokens) < 2) {
        range_pairs = (min(q_grid), max(q_grid))
    }
    else {
        real scalar n_pairs
        n_pairs = floor(cols(range_tokens) / 2)
        range_pairs = J(n_pairs, 2, .)
        real scalar idx
        idx = 1
        for (real scalar r = 1; r <= n_pairs; r++) {
            range_pairs[r,1] = range_tokens[idx]
            range_pairs[r,2] = range_tokens[idx+1]
            idx = idx + 2
        }
    }

    // For simplicity, use the first range if multiple provided
    real scalar range_lo, range_hi
    range_lo = range_pairs[1,1]
    range_hi = range_pairs[1,2]

    real colvector range_idx
    range_idx = J(0, 1, .)
    for (real scalar j = 1; j <= nq; j++) {
        if (q_grid[j] >= range_lo & q_grid[j] <= range_hi) {
            range_idx = range_idx \\ j
        }
    }

    if (rows(range_idx) == 0) {
        range_idx = (1..nq)'
    }

    real matrix pvals
    pvals = J(1, 3, .)

    if (do_nullity) {
        real scalar test_stat_null
        test_stat_null = max(abs(tau_orig[range_idx]'))

        real colvector supvals_null
        supvals_null = J(B, 1, 0)
        for (real scalar b = 1; b <= B; b++) {
            supvals_null[b] = max(abs(tau_boot[b, range_idx]'))
        }

        real colvector sup_null_sorted
        sup_null_sorted = sort(supvals_null, 1)
        real scalar idx_null
        idx_null = ceil(prob * (B + 1))
        if (idx_null < 1) idx_null = 1
        if (idx_null > B) idx_null = B
        real scalar crit_null
        crit_null = sup_null_sorted[idx_null]

        real scalar p_null
        p_null = mean(supvals_null :>= test_stat_null)
        pvals[1,1] = p_null
    }

    if (do_homogeneity) {
        real scalar test_stat_homo
        real rowvector range_tau
        range_tau = tau_orig[range_idx]
        real scalar mbar
        mbar = mean(range_tau)
        test_stat_homo = max(abs(range_tau :- mbar))

        real colvector supvals_homo
        supvals_homo = J(B, 1, 0)
        for (real scalar b = 1; b <= B; b++) {
            real rowvector draw
            draw = tau_boot[b, range_idx]
            real scalar draw_mean
            draw_mean = mean(draw)
            supvals_homo[b] = max(abs(draw :- draw_mean))
        }

        real colvector sup_homo_sorted
        sup_homo_sorted = sort(supvals_homo, 1)
        real scalar idx_homo
        idx_homo = ceil(prob * (B + 1))
        if (idx_homo < 1) idx_homo = 1
        if (idx_homo > B) idx_homo = B

        real scalar crit_homo
        crit_homo = sup_homo_sorted[idx_homo]

        real scalar p_homo
        p_homo = mean(supvals_homo :>= test_stat_homo)
        pvals[1,2] = p_homo
    }

    if (do_gini) {
        real scalar gini_above, gini_below, gini_diff
        gini_above = r3d_gini_from_quantile(q_grid, int_plus)
        gini_below = r3d_gini_from_quantile(q_grid, int_minus)
        gini_diff = gini_above - gini_below

        real colvector gini_boot
        gini_boot = J(B, 1, 0)
        for (real scalar b = 1; b <= B; b++) {
            real rowvector boot_plus, boot_minus
            boot_plus = int_plus + nu_plus_store[b,]
            boot_minus = int_minus + nu_minus_store[b,]
            gini_boot[b] = r3d_gini_from_quantile(q_grid, boot_plus) - ///
                           r3d_gini_from_quantile(q_grid, boot_minus)
        }

        real scalar test_stat_gini
        test_stat_gini = abs(gini_diff)

        real colvector gini_sorted
        gini_sorted = sort(abs(gini_boot), 1)
        real scalar idx_gini
        idx_gini = ceil(prob * (B + 1))
        if (idx_gini < 1) idx_gini = 1
        if (idx_gini > B) idx_gini = B

        real scalar crit_gini
        crit_gini = gini_sorted[idx_gini]

        real scalar p_gini
        p_gini = mean(abs(gini_boot) :>= test_stat_gini)
        pvals[1,3] = p_gini
    }

    st_matrix(bs_results_name, tau_boot)
    st_matrix(pvals_name, pvals)
    st_matrix("r(cb_lower)", cb_lower)
    st_matrix("r(cb_upper)", cb_upper)
}

// ============================================================================
// GINI COEFFICIENT
// ============================================================================

// Calculate Gini coefficient from quantile function
real scalar r3d_gini_from_quantile(real rowvector quantiles, real rowvector qfunc)
{
    real scalar nq, i, mean_val, area, gini

    nq = cols(quantiles)
    if (nq < 2) return(.)

    // Approximate mean from quantile function
    mean_val = 0
    for (i = 1; i < nq; i++) {
        mean_val = mean_val + 0.5 * (qfunc[i+1] + qfunc[i]) * (quantiles[i+1] - quantiles[i])
    }

    if (mean_val <= 0) return(0)

    // Calculate Lorenz curve
    real rowvector lorenz
    lorenz = J(1, nq, 0)
    real scalar cumulative
    cumulative = 0

    for (i = 1; i < nq; i++) {
        cumulative = cumulative + 0.5 * (qfunc[i+1] + qfunc[i]) * (quantiles[i+1] - quantiles[i])
        lorenz[i+1] = cumulative / mean_val
    }

    // Area under Lorenz curve
    area = 0
    for (i = 1; i < nq; i++) {
        area = area + 0.5 * (lorenz[i] + lorenz[i+1]) * (quantiles[i+1] - quantiles[i])
    }

    // Gini coefficient
    gini = 1 - 2 * area

    return(max((0, min((1, gini)))))
}

// Test for Gini coefficient difference
void r3d_test_gini(string scalar xvar, string scalar Qmat,
                   real scalar cutoff, string scalar touse,
                   string scalar pval_name)
{
    real matrix X, Q
    real scalar n, nq, gini_plus, gini_minus, diff
    real rowvector quantiles, qfunc_plus, qfunc_minus

    // Get data
    st_view(X, ., xvar, touse)
    Q = st_matrix(Qmat)
    n = rows(X)
    nq = cols(Q)

    // Create quantile grid
    quantiles = J(1, nq, 0)
    real scalar i
    for (i = 1; i <= nq; i++) {
        quantiles[i] = i / (nq + 1)
    }

    // Compute average quantile functions above and below cutoff
    real colvector idx_plus, idx_minus
    idx_plus = (X :>= 0)
    idx_minus = (X :< 0)

    qfunc_plus = J(1, nq, 0)
    qfunc_minus = J(1, nq, 0)

    real scalar j
    for (j = 1; j <= nq; j++) {
        qfunc_plus[j] = mean(select(Q[,j], idx_plus))
        qfunc_minus[j] = mean(select(Q[,j], idx_minus))
    }

    // Calculate Gini coefficients
    gini_plus = r3d_gini_from_quantile(quantiles, qfunc_plus)
    gini_minus = r3d_gini_from_quantile(quantiles, qfunc_minus)

    // Difference
    diff = gini_plus - gini_minus

    // Store results (p-value would need bootstrap for proper inference)
    st_numscalar(pval_name, 0.5)  // Placeholder
    st_numscalar("r(gini_plus)", gini_plus)
    st_numscalar("r(gini_minus)", gini_minus)
    st_numscalar("r(gini_diff)", diff)
}

end
