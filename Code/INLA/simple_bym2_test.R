#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Minimal, self-contained INLA BYM2 model test for HPC validation.
# This script does not require any external data files.

cat("=====================================================\n")
cat("=      Minimal INLA BYM2 HPC Validation Test      =\n")
cat("=====================================================\n\n")

# 1. Load INLA library
cat("[1/5] Loading INLA library...\n")
suppressMessages(library(INLA))
cat("      ✓ INLA loaded successfully.\n\n")

# 2. Create a spatial structure (10x10 grid = 100 areas)
cat("[2/5] Creating a 10x10 virtual spatial grid (100 areas)...\n")
n_side <- 10
n_areas <- n_side * n_side
adj_matrix <- matrix(0, nrow = n_areas, ncol = n_areas)

for (i in 1:n_areas) {
  # Neighbor to the right
  if (i %% n_side != 0) {
    adj_matrix[i, i + 1] <- 1
    adj_matrix[i + 1, i] <- 1
  }
  # Neighbor below
  if (i <= (n_areas - n_side)) {
    adj_matrix[i, i + n_side] <- 1
    adj_matrix[i + n_side, i] <- 1
  }
}
# Create INLA graph object from the adjacency matrix
g <- inla.read.graph(adj_matrix)
cat("      ✓ Spatial graph created successfully.\n\n")

# 3. Generate synthetic data for the 100 areas
cat("[3/5] Generating synthetic data for the model...\n")
set.seed(12345) # for reproducibility
alpha <- -1.0   # True intercept
beta <- 0.5     # True effect for the covariate
x <- rnorm(n_areas) # A random covariate
E <- runif(n_areas, min = 10, max = 50) # Expected counts (e.g., from population)

# Linear predictor and Poisson response
eta <- alpha + beta * x
lambda <- E * exp(eta)
y <- rpois(n_areas, lambda)

# Create the data frame for INLA
df <- data.frame(
  id = 1:n_areas,
  y = y,
  x = x,
  E = E
)
cat("      ✓ Synthetic data generated.\n\n")

# 4. Define and fit the BYM2 model
cat("[4/5] Defining and fitting the INLA BYM2 model...\n")
formula <- y ~ 1 + x + f(id, model = "bym2", graph = g)

# Set INLA to use a reasonable number of threads
inla.setOption(num.threads = "4:1")

# Run the model
result <- tryCatch({
  inla(formula,
       family = "poisson",
       data = df,
       E = E, # Use E as the exposure/offset term
       control.predictor = list(compute = TRUE))
}, error = function(e) {
  cat("      ❌ ERROR during model fitting:\n")
  print(e)
  return(NULL)
})
cat("      ✓ Model fitting process completed.\n\n")

# 5. Print the results summary
cat("[5/5] Displaying model summary...\n")
if (!is.null(result)) {
  print(summary(result))
  cat("\n=====================================================\n")
  cat("=  ✅ TEST COMPLETED SUCCESSFULLY!                =\n")
  cat("=====================================================\n")
} else {
  cat("\n=====================================================\n")
  cat("=  ❌ TEST FAILED!                                =\n")
  cat("=====================================================\n")
}
