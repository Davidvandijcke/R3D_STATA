#!/usr/bin/env Rscript
#
# generate_reference_data.R
# Generates reference datasets and R3D outputs for Stata equivalence testing.
# Run from tests/ directory: Rscript generate_reference_data.R
#

library(R3D)

cat("=== Generating R3D reference data for Stata equivalence tests ===\n")

# ============================================================================
# Helper: create data and save to CSV for Stata
# ============================================================================
save_data_for_stata <- function(X, Y_list, q_grid, filename, T_vec = NULL) {
  n <- length(X)
  nq <- length(q_grid)

  # Compute quantile matrix
  Qmat <- R3D:::.compute_empirical_qmat(Y_list, q_grid)

  # Build data frame: X, T (if fuzzy), Q1..Qnq
  df <- data.frame(X = X)
  if (!is.null(T_vec)) {
    df$T_treat <- T_vec
  }
  for (j in seq_len(nq)) {
    df[[paste0("Q", j)]] <- Qmat[, j]
  }

  write.csv(df, filename, row.names = FALSE)
  cat("  Saved:", filename, "(", n, "obs x", ncol(df), "vars)\n")
  return(Qmat)
}

save_results <- function(fit, filename) {
  nq <- length(fit$q_grid)
  df <- data.frame(
    q = fit$q_grid,
    tau = as.numeric(fit$tau),
    bw_num = as.numeric(fit$bandwidths$h_star_num),
    int_plus = as.numeric(fit$int_plus),
    int_minus = as.numeric(fit$int_minus)
  )
  if (!is.null(fit$bandwidths$h_star_den) && length(fit$bandwidths$h_star_den) > 0) {
    df$bw_den <- rep(fit$bandwidths$h_star_den, nq)
  }

  # Add bootstrap results if available
  if (!is.null(fit$boot_out)) {
    df$cb_lower <- as.numeric(fit$boot_out$cb_lower)
    df$cb_upper <- as.numeric(fit$boot_out$cb_upper)
  }

  write.csv(df, filename, row.names = FALSE)
  cat("  Saved:", filename, "\n")
}

# ============================================================================
# Case 1: Sharp-Simple, 20 quantiles, Epanechnikov kernel
# ============================================================================
cat("\n--- Case 1: Sharp-Simple, nq=20, Epanechnikov ---\n")
set.seed(42)
n1 <- 200
X1 <- runif(n1, -1, 1)
Y_list1 <- lapply(seq_len(n1), function(i) {
  mu <- 5 + 0.5 * X1[i] + 3 * (X1[i] >= 0)
  rnorm(50, mean = mu, sd = 0.5)
})
q_grid_20 <- (1:20) / 21

Qmat1 <- save_data_for_stata(X1, Y_list1, q_grid_20, "ref_data_sharp_simple_20.csv")

fit1 <- r3d(X1, Y_list1, cutoff = 0, method = "simple", p = 2,
            q_grid = q_grid_20, kernel_fun = "epanechnikov")
save_results(fit1, "ref_results_sharp_simple_20.csv")

# ============================================================================
# Case 2: Sharp-Simple, 99 quantiles, Triangular kernel
# ============================================================================
cat("\n--- Case 2: Sharp-Simple, nq=99, Triangular ---\n")
set.seed(42)
q_grid_99 <- (1:99) / 100

Qmat2 <- save_data_for_stata(X1, Y_list1, q_grid_99, "ref_data_sharp_simple_99.csv")

fit2 <- r3d(X1, Y_list1, cutoff = 0, method = "simple", p = 2,
            q_grid = q_grid_99, kernel_fun = "triangular")
save_results(fit2, "ref_results_sharp_simple_99.csv")

# ============================================================================
# Case 3: Sharp-Frechet, 20 quantiles, Epanechnikov kernel
# ============================================================================
cat("\n--- Case 3: Sharp-Frechet, nq=20, Epanechnikov ---\n")
set.seed(42)
fit3 <- r3d(X1, Y_list1, cutoff = 0, method = "frechet", p = 2,
            q_grid = q_grid_20, kernel_fun = "epanechnikov")
save_results(fit3, "ref_results_sharp_frechet_20.csv")

# ============================================================================
# Case 4: Sharp-Simple, Uniform kernel
# ============================================================================
cat("\n--- Case 4: Sharp-Simple, nq=20, Uniform ---\n")
set.seed(42)
fit4 <- r3d(X1, Y_list1, cutoff = 0, method = "simple", p = 2,
            q_grid = q_grid_20, kernel_fun = "uniform")
save_results(fit4, "ref_results_sharp_simple_uniform.csv")

# ============================================================================
# Case 5: Fuzzy-Simple, 20 quantiles
# ============================================================================
cat("\n--- Case 5: Fuzzy-Simple, nq=20, Epanechnikov ---\n")
set.seed(42)
n5 <- 200
X5 <- runif(n5, -1, 1)
p_treat <- 0.1 + 0.8 * (X5 >= 0)
T5 <- rbinom(n5, 1, p_treat)
Y_list5 <- lapply(seq_len(n5), function(i) {
  mu <- 5 + 0.5 * X5[i] + 3 * T5[i]
  rnorm(50, mean = mu, sd = 0.5)
})

Qmat5 <- save_data_for_stata(X5, Y_list5, q_grid_20, "ref_data_fuzzy_simple_20.csv", T_vec = T5)

fit5 <- r3d(X5, Y_list5, T = T5, cutoff = 0, method = "simple", p = 2,
            q_grid = q_grid_20, kernel_fun = "epanechnikov", fuzzy = TRUE)
save_results(fit5, "ref_results_fuzzy_simple_20.csv")

# ============================================================================
# Case 6: Fuzzy-Frechet, 20 quantiles
# ============================================================================
cat("\n--- Case 6: Fuzzy-Frechet, nq=20, Epanechnikov ---\n")
set.seed(42)
fit6 <- r3d(X5, Y_list5, T = T5, cutoff = 0, method = "frechet", p = 2,
            q_grid = q_grid_20, kernel_fun = "epanechnikov", fuzzy = TRUE)
save_results(fit6, "ref_results_fuzzy_frechet_20.csv")

# ============================================================================
# Case 7: Bootstrap with pre-generated multipliers (sharp-simple)
# ============================================================================
cat("\n--- Case 7: Bootstrap, sharp-simple, nq=20 ---\n")
set.seed(42)
fit7 <- r3d(X1, Y_list1, cutoff = 0, method = "simple", p = 2,
            q_grid = q_grid_20, kernel_fun = "epanechnikov",
            boot = TRUE, boot_reps = 500, alpha = 0.05,
            test = c("nullity", "homogeneity"))
save_results(fit7, "ref_results_bootstrap_sharp.csv")

# Save p-values from test_results
if (!is.null(fit7$boot_out$test_results)) {
  test_names <- c()
  p_vals <- c()
  for (tname in names(fit7$boot_out$test_results)) {
    for (rname in names(fit7$boot_out$test_results[[tname]])) {
      tr <- fit7$boot_out$test_results[[tname]][[rname]]
      test_names <- c(test_names, paste0(tname, ":", rname))
      p_vals <- c(p_vals, tr$p_value)
    }
  }
  pv <- data.frame(test = test_names, p_value = p_vals)
  write.csv(pv, "ref_pvalues_bootstrap_sharp.csv", row.names = FALSE)
  cat("  Saved: ref_pvalues_bootstrap_sharp.csv\n")
}

# ============================================================================
# Save density estimate for verification
# ============================================================================
cat("\n--- Density estimate verification ---\n")
Xc <- X1 - 0
sigma_X <- sd(Xc)
h_bw <- 1.06 * sigma_X * n1^(-1/5)
kernel_tri <- function(u) pmax(0, 1 - abs(u))
f_X_hat <- mean(kernel_tri(Xc / h_bw)) / h_bw
cat("  f_X_hat (triangular):", f_X_hat, "\n")
cat("  h_bw:", h_bw, "\n")
cat("  sigma_X:", sigma_X, "\n")

density_df <- data.frame(
  f_X_hat_tri = f_X_hat,
  h_bw = h_bw,
  sigma_X = sigma_X
)
write.csv(density_df, "ref_density.csv", row.names = FALSE)

cat("\n=== Reference data generation complete ===\n")
