dynamic_sampler <- nimbleFunction(
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    # Parâmetros/objetos passados via control
    E         <- control$E
    x         <- control$x
    a0        <- control$a0
    b0        <- control$b0
    w         <- control$w
    n_regions <- control$n_regions
    n_times   <- control$n_times
  },
  
  run = function() {
    ## Matrizes auxiliares locais
    g_b       <- matrix(0, n_regions, n_times)
    sum_y     <- matrix(0, n_regions, n_times)
    sum_g_b   <- matrix(0, n_regions, n_times)
    a         <- matrix(0, n_regions, n_times)
    b         <- matrix(0, n_regions, n_times)
    
    ## --------- ESCALARES (tipo double) para evitar NimArr<1,double> ----------
    declare(E_rt,        double(0))
    declare(eps_rt,      double(0))
    declare(x_rt,        double(0))
    declare(beta_val,    double(0))
    declare(g_tmp,       double(0))
    declare(sum_g_tmp,   double(0))
    declare(b_tmp,       double(0))
    declare(y_rt,        double(0))
    declare(a_local,     double(0))
    declare(b_local,     double(0))
    declare(lambda_aux,  double(0))
    declare(t_idx,       integer(0))
    ## -------------------------------------------------------------------------
    
    ## t = 1
    for (r in 1:n_regions) {
      a[r, 1] <- a0
      b[r, 1] <- b0
      
      E_rt     <- E[r, 1]
      eps_rt   <- model$epsilon[r, 1]
      x_rt     <- x[r, 1]
      beta_val <- model$beta[1]
      
      g_tmp        <- E_rt * eps_rt * exp(beta_val * x_rt)
      g_b[r, 1]    <- g_tmp
      sum_g_b[r,1] <- g_tmp
      
      y_rt         <- model$Y[r, 1]
      sum_y[r, 1]  <- y_rt
      
      model$lambda[r, 1] <<- rgamma(1, shape = a[r, 1], scale = 1 / b[r, 1])
    }
    
    ## Forward: t >= 2
    for (t in 2:n_times) {
      for (r in 1:n_regions) {
        y_rt     <- model$Y[r, t]
        E_rt     <- E[r, t]
        eps_rt   <- model$epsilon[r, t]
        x_rt     <- x[r, t]
        beta_val <- model$beta[1]
        
        sum_y[r, t] <- sum_y[r, t - 1] + y_rt
        a[r, t]     <- w * a[r, t - 1] + sum_y[r, t]
        
        g_tmp        <- E_rt * eps_rt * exp(beta_val * x_rt)
        sum_g_tmp    <- sum_g_b[r, t - 1] + g_tmp
        b_tmp        <- w * b[r, t - 1] + sum_g_tmp
        
        g_b[r, t]     <- g_tmp
        sum_g_b[r, t] <- sum_g_tmp
        b[r, t]       <- b_tmp
        
        model$lambda[r, t] <<- rgamma(1, shape = a[r, t], scale = 1 / b[r, t])
      }
    }
    
    ## Backward smoothing
    for (r in 1:n_regions) {
      lambda_aux <- model$lambda[r, n_times]  # só reutilizando a variável
      # inicializa o último tempo
      # usamos g_b como temporário só para garantir shape; mas melhor criar matriz própria:
    }
    lambda_smooth <- matrix(0, n_regions, n_times)
    for (r in 1:n_regions) {
      lambda_smooth[r, n_times] <- model$lambda[r, n_times]
    }
    
    for (j in 1:(n_times - 1)) {
      t_idx <- n_times - j
      for (r in 1:n_regions) {
        a_local    <- a[r, t_idx]
        b_local    <- b[r, t_idx]
        lambda_aux <- rgamma(1, shape = (1 - w) * a_local, scale = 1 / b_local)
        lambda_smooth[r, t_idx] <- lambda_aux + w * lambda_smooth[r, t_idx + 1]
      }
    }
    
    ## Atualização elemento-a-elemento no modelo
    for (r in 1:n_regions) {
      for (t in 1:n_times) {
        model$lambda[r, t] <<- lambda_smooth[r, t]
      }
    }
    
    copy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
  },
  
  methods = list(reset = function() {})
)
