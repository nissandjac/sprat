
par <- sas_2024$reps$par.fixed

alpha_seq <- seq(par["pin"] - 2, par["pin"] + 2, length.out = 40)
beta_seq  <- seq(par["logQ"][1]  - 2, par["logQ"][1]  + 2, length.out = 40)

grid <- expand.grid(alpha = alpha_seq, beta = beta_seq)
nll  <- numeric(nrow(grid))

for (i in seq_len(nrow(grid))) {
  par_test <- sas_2024$reps$par.fixed
  par_test["pin"] <- grid$alpha[i]
  par_test["logQ"][1]  <- grid$beta[i]
  nll[i] <- sas_2024$obj$fn(par_test)  # NO re-optimisation: pure surface
}

mat_nll <- matrix(nll, nrow = length(alpha_seq), byrow = FALSE)
contour(alpha_seq, beta_seq, mat_nll,
        xlab = "Power law", ylab = "Q age 0", main = "NLL surface")


