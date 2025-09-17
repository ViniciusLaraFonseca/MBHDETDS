dynamic_sampler <- nimbleFunction(
  contains = sampler_BASE,  # Inherit from sampler_BASE
  
  setup = function(model, mvSaved, target, control) {
    # Access model parameters and data, storing them in the scope of the nimbleFunction
    epsilon <<- model$epsilon  # Effect of covariates on the model
    E <- modelConstants$E              # Constant values for each region and time
    beta <<- model$beta        # Regression coefficients
    y <<- model$Y              # Observed counts
    a0 <- modelConstants$a0            # Initial value for parameter a
    b0 <- modelConstants$b0            # Initial value for parameter b
    w <- modelConstants$w              # Weighting factor
    n_regions <- modelConstants$n_regions  # Number of regions
    n_times <- modelConstants$n_times      # Number of time points
  },
  
  run = function() {
    # Initialize matrices for parameters a, b, g_b, and lambda
    a <- matrix(0, nrow = n_regions, ncol = n_times)  # Matrix for parameter a
    b <- matrix(0, nrow = n_regions, ncol = n_times)  # Matrix for parameter b
    g_b <- matrix(0, nrow = n_regions, ncol = n_times)  # Matrix for auxiliary parameter g_b
    lambda <- matrix(0, nrow = n_regions, ncol = n_times)  # Matrix for lambda
    
    # Set initial values for the first time point
    a[, 1] <- a0  # Set initial a values
    b[, 1] <- b0  # Set initial b values
    lambda[, 1] <- rgamma(n_regions, shape = a[, 1], scale = 1 / b[, 1])  # Initialize lambda with gamma distribution
    
    # Forward updates for each time point
    for (t in 2:n_times) {  # Loop over time points starting from the second
      for (r in 1:n_regions) {  # Loop over each region
        sum_y <- sum(y[r, 1:t])  # Calculate the cumulative sum of observed values for region r up to time t
        a[r, t] <- w * a[r, t - 1] + sum_y  # Update parameter a using the previous value and cumulative sum
      }
      
      for (r in 1:n_regions) {  # Loop over each region for g_b calculation
        g_b[r, t] <- E[r, t] * epsilon[r, t] * exp(beta %*% t(x[, t, ]))  # Calculate g_b for region r at time t
      }
      
      for (r in 1:n_regions) {  # Loop over each region to update b
        sum_g_b <- sum(g_b[r, 1:t])  # Calculate the cumulative sum of g_b for region r up to time t
        b[r, t] <- w * b[r, t - 1] + sum_g_b  # Update parameter b using the previous value and cumulative sum
      }
      
      lambda[, t] <- rgamma(n_regions, shape = a[, t], scale = 1 / b[, t])  # Draw new lambda values from gamma distribution
    }
    
    #Initialize matrices for auxiliary and smoothed lambda values
    lambda_aux <- matrix(0, nrow = n_regions, ncol = n_times)  # Auxiliary lambda values
    lambda_smooth <- matrix(0, nrow = n_regions, ncol = n_times)  # Smoothed lambda values
    
    #Set the last value of smoothed lambda to the last value of lambda
    lambda_smooth[, n_times] <- lambda[, n_times]  # Initialize the last column of smoothed lambda
    
    #Smoothing lambda values from the last time point to the first
    for (j in 1:(n_times - 1)) {  # Loop over time points in reverse
      for (r in 1:n_regions) {  # Loop over each region
        #Generate auxiliary lambda values
        lambda_aux[r, n_times - j] <- rgamma(1, shape = (1 - w) * a[r, n_times - j], scale = 1 / b[r, n_times - j])  # Draw from gamma distribution
        #Smooth the lambda values
        lambda_smooth[r, n_times - j] <- lambda_aux[r, n_times - j] + w * lambda[r, n_times - j + 1]  # Combine auxiliary and previous lambda
      }
    }
    
    #Update the lambda values in the model with the smoothed lambda
    for (r in 1:n_regions) {  # Loop over each region
      model$lambda[r, ] <<- lambda_smooth[r, ]  # Store smoothed lambda back in the model
    }
    
    copy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)  # Save the updated model values
  }
)
