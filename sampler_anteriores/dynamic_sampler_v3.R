if(!require(nimble)){ install.packages("nimble"); require(nimble)}
dynamic_sampler <- nimbleFunction(
  contains = sampler_BASE,  # Inherit from sampler_BASE
  
  setup = function(model, mvSaved, target, control) {
    # Access model parameters and data, storing them in the scope of the nimbleFunction
    epsilon <- model$epsilon   # Effect of covariates on the model
    E <-control$E      # Constant values for each region and time
    beta <- model$beta         # Regression coefficients
    Y <- model$Y               # Observed counts
    a0 <-control$a0      # Initial value for parameter a
    b0 <- control$b0      # Initial value for parameter b
    w <- control$w        # Weighting factor
    x <- control$x
    n_regions <- control$n_regions  # Number of regions
    n_times <- control$n_times      # Number of time points
  },
  
  run = function() {
    ignore <- model$calculate()
    # Initialize matrices LOCALLY using <-
    g_b <- matrix(0, n_regions, n_times)
    lambda <- matrix(0, n_regions, n_times)
    sum_y <- matrix(0, n_regions, n_times)
    sum_g_b <- matrix(0, n_regions, n_times)
    
    # Initialize time 1 properly
    
    for (r in 1:n_regions) {
      model$a[r, 1] <<- a0
      model$b[r, 1] <<- b0
      g_b[r, 1] <- E[r, 1] * model$epsilon[r, 1] * exp(inprod(model$beta , x[r, 1, ]))
      sum_g_b[r, 1] <- g_b[r, 1]
      sum_y[r, 1] <- model$Y[r, 1]
      lambda[r, 1] <- rgamma(1, shape = model$a[r, 1], scale = 1 / model$b[r, 1])
    }
    
    # Forward loop (t >= 2)
    for (t in 2:n_times) {
      for (r in 1:n_regions) {
        sum_y[r, t] <- sum_y[r, t-1] + model$Y[r, t]  # Recursive cumulative sum
        model$a[r, t] <<- w * model$a[r, t-1] + sum_y[r, t]  # Use <- for assignment
      }
      for (r in 1:n_regions) {
        g_b[r, t] <- E[r, t] * model$epsilon[r, t] * exp(inprod(model$beta , x[r, t, ]))
        sum_g_b[r, t] <- sum_g_b[r, t-1] + g_b[r, t]  # Recursive cumulative
        model$b[r, t] <<- w * model$b[r, t-1] + sum_g_b[r, t]  # Use <- for assignment
      }
      for (r in 1:n_regions) {
        model$lambda[r, t] <<- rgamma(1, shape = model$a[r, t], scale = 1 / model$b[r, t])  # rgamma not ~
      }
    }
    
    # Smoothing (unchanged but use <-)
    lambda_aux <- matrix(0, n_regions, n_times)
    lambda_smooth <- matrix(0, n_regions, n_times)
    lambda_smooth[, n_times] <- model$lambda[, n_times]
    for (j in 1:(n_times - 1)) {
      for (r in 1:n_regions) {
        lambda_aux[r, n_times - j] <- rgamma(1, shape = (1 - w) * model$a[r, n_times - j], scale = 1 / model$b[r, n_times - j])
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
