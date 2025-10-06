#!/usr/bin/env Rscript
library(brms)

# Create simple test data
set.seed(123)
df <- data.frame(
  y = rnorm(100, 50, 10),
  x = rnorm(100),
  group = factor(rep(1:10, 10))
)

# Test simple prior specification
priors <- c(
  set_prior("normal(50, 10)", class = "Intercept"),
  set_prior("normal(0, 2)", class = "b"),
  set_prior("student_t(3, 0, 10)", class = "sigma")
)

# Try to fit simple model
tryCatch({
  fit <- brm(y ~ x, data = df, prior = priors, chains = 1, iter = 100, cores = 1, refresh = 0)
  cat("Simple brms model works!\n")
}, error = function(e) {
  cat("Error in brms model:", e$message, "\n")
})