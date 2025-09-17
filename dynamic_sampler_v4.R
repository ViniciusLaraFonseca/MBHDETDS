dynamic_sampler <- nimbleFunction(
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    beta <- model$beta         # Regression coefficients
    Y <- model$Y  
    epsilon <- model$epsilon
    E <- control$E
    w <- control$w
    x <- control$x
    a0 <- control$a0
    b0 <- control$b0
    n_regions <- control$n_regions
    n_times <- control$n_times
  },
  
  run = function() {
    # Matrizes auxiliares locais
    g_b <- matrix(0, n_regions, n_times)
    sum_y <- matrix(0, n_regions, n_times)
    sum_g_b <- matrix(0, n_regions, n_times)
    a <- matrix(0, n_regions, n_times)
    b <- matrix(0, n_regions, n_times)
    
    # Inicialização em t = 1
    for (r in 1:n_regions) {
      a[r, 1] <- a0
      b[r, 1] <- b0
      g_b[r, 1] <- E[r, 1] * model$epsilon[r, 1] * exp(model$beta * x[r,1])
      sum_g_b[r, 1] <- g_b[r, 1]
      sum_y[r, 1] <- model$Y[r, 1]
      model$lambda[r, 1] <<- rgamma(1, shape = a[r, 1], scale = 1 / b[r, 1])
    }
    
    # Forward (t >= 2)
    for (t in 2:n_times) {
      for (r in 1:n_regions) {
        sum_y[r, t] <- sum_y[r, t - 1] + model$Y[r, t]
        a[r, t] <- w * a[r, t - 1] + sum_y[r, t]
        
        g_b[r, t] <- E[r, t] * model$epsilon[r, t] * exp(model$beta * x[r,t])
        sum_g_b[r, t] <- sum_g_b[r, t - 1] + g_b[r, t]
        b[r, t] <- w * b[r, t - 1] + sum_g_b[r, t]
        
        model$lambda[r, t] <<- rgamma(1, shape = a[r, t], scale = 1 / b[r, t])
      }
    }
    
    # Backward smoothing
    lambda_smooth <- matrix(0, n_regions, n_times)
    lambda_smooth[, n_times] <- model$lambda[, n_times]
    
    for (j in 1:(n_times - 1)) {
      for (r in 1:n_regions) {
        lambda_aux <- rgamma(1, shape = (1 - w) * a[r, n_times - j], scale = 1 / b[r, n_times - j])
        lambda_smooth[r, n_times - j] <- lambda_aux + w * lambda_smooth[r, n_times - j + 1]
      }
    }
    
    # Atualizar no modelo
    for (r in 1:n_regions) {
      model$lambda[r, ] <<- lambda_smooth[r, ]
    }
    
    copy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
  },
  
  methods = list(reset = function() {})
)
