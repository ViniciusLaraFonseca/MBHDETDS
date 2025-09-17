if(!require(nimble)){ install.packages("nimble"); require(nimble)}
dynamic_sampler <- nimbleFunction(
  contains = sampler_BASE,  # Inherit from sampler_BASE
  
  setup = function(model, mvSaved, target, control) {
    # Access model parameters and data, storing them in the scope of the nimbleFunction
    epsilon <- model$epsilon   # Effect of covariates on the model
    E <- modelConstants$E      # Constant values for each region and time
    beta_values <- model$beta         # Regression coefficients
    Y <- model$Y               # Observed counts
    a0 <- modelConstants$a0      # Initial value for parameter a
    b0 <- modelConstants$b0      # Initial value for parameter b
    w <- modelConstants$w        # Weighting factor
    x <- modelConstants$x
    n_regions <- modelConstants$n_regions  # Number of regions
    n_times <- modelConstants$n_times      # Number of time points
  },
  
  run = function() {
    cat("Running dynamic sampler for lambda...\n")
    # Initialize matrices LOCALLY using <-
    a <- matrix(0, n_regions, n_times)
    b <- matrix(0, n_regions, n_times)
    g_b <- matrix(0, n_regions, n_times)
    lambda <- matrix(0, n_regions, n_times)
    sum_y <- matrix(0, n_regions, n_times)
    sum_g_b <- matrix(0, n_regions, n_times)
    
    # Initialize time 1 properly
   
    for (r in 1:n_regions) {
      a[r, 1] <- a0
      b[r, 1] <- b0
      g_b[r, 1] <- E[r, 1] * epsilon[r, 1] * exp(inprod(beta_values, x[r, 1, ]))
      sum_g_b[r, 1] <- g_b[r, 1]
      sum_y[r, 1] <- Y[r, 1]
      lambda[r, 1] <- rgamma(1, shape = a[r, 1], scale = 1 / b[r, 1])
    }
    
    # Forward loop (t >= 2)
    for (t in 2:n_times) {
      for (r in 1:n_regions) {
        sum_y[r, t] <- sum_y[r, t-1] + Y[r, t]  # Recursive cumulative sum
        model$a[r, t] <<- w * a[r, t-1] + sum_y[r, t]  # Use <- for assignment
      }
      for (r in 1:n_regions) {
        g_b[r, t] <- E[r, t] * epsilon[r, t] * exp(inprod(beta_values, x[r, t, ]))
        sum_g_b[r, t] <- sum_g_b[r, t-1] + g_b[r, t]  # Recursive cumulative
        model$b[r, t] <<- w * b[r, t-1] + sum_g_b[r, t]  # Use <- for assignment
      }
      for (r in 1:n_regions) {
        lambda[r, t] <- rgamma(1, shape = a[r, t], scale = 1 / b[r, t])  # rgamma not ~
      }
    }
    
    # Smoothing (unchanged but use <-)
    lambda_aux <- matrix(0, n_regions, n_times)
    lambda_smooth <- matrix(0, n_regions, n_times)
    lambda_smooth[, n_times] <- lambda[, n_times]
    for (j in 1:(n_times - 1)) {
      for (r in 1:n_regions) {
        lambda_aux[r, n_times - j] <- rgamma(1, shape = (1 - w) * a[r, n_times - j], scale = 1 / b[r, n_times - j])
        lambda_smooth[r, n_times - j] <- lambda_aux[r, n_times - j] + w * lambda_smooth[r, n_times - j + 1]
      }
    }
    
    # Update model (only place <<- is needed)
    for (r in 1:n_regions) {
      model$lambda[r, ] <<- lambda_smooth[r, ]
    }
    copy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
  },
  methods = list(reset = function() {})
)
