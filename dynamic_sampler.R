dynamic_sampler <- nimbleFunction(
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    # Extrair dimensões e parâmetros
    dims_Y <- dim(model$Y)
    n_regions <- dims_Y[1]
    n_times   <- dims_Y[2]
    p         <- dim(model$x)[3]
    w         <- control$w
    a0        <- control$a0
    b0        <- control$b0
    
    # Buffers para armazenar os parâmetros da filtragem
    at_buf  <- nimMatrix(nrow = n_regions, ncol = n_times + 1, init = 0, type = 'double')
    bt_buf  <- nimMatrix(nrow = n_regions, ncol = n_times + 1, init = 0, type = 'double')
    
    # Nós a serem calculados após a amostragem
    calcNodes   <- model$getDependencies(target, self = FALSE)
    targetNodes <- model$expandNodeNames(target)
  },
  
  run = function() {
    # --- Declaração explícita de tipos ---
    declare(i, integer())
    declare(t, integer())
    declare(k, integer())
    declare(prod_val, double())
    declare(att_t, double())
    declare(btt_t, double())
    declare(shape_tmp, double())
    declare(rate_tmp, double())
    declare(lambda_futuro, double())
    declare(nu, double())
    
    # --- Forward Filtering ---
    for(i in 1:n_regions) {
      at_buf[i, 1] <<- a0
      bt_buf[i, 1] <<- b0
      
      for(t in 1:n_times) {
        att_t <- w * at_buf[i, t]
        btt_t <- w * bt_buf[i, t]
        
        at_buf[i, t+1] <<- att_t + model$Y[i, t]
        
        prod_val <- 0
        for(k in 1:p) {
          prod_val <- prod_val + model$x[i, t, k] * model$beta[k]
        }
        
        bt_buf[i, t+1] <<- btt_t + model$E[i, t] * model$epsilon[i] * exp(prod_val)
      }
    }
    
    # --- Backward Sampling ---
    for(i in 1:n_regions) {
      shape_tmp <- at_buf[i, n_times + 1]
      rate_tmp  <- bt_buf[i, n_times + 1]
      model$lambda[i, n_times] <<- rgamma(1, shape = shape_tmp, rate = rate_tmp)
      
      for(t in n_times:2) {
        lambda_futuro <- model$lambda[i, t]
        
        shape_tmp <- (1 - w) * at_buf[i, t] 
        rate_tmp  <- bt_buf[i, t]
        
        nu <- rgamma(1, shape = shape_tmp, rate = rate_tmp)
        
        model$lambda[i, t-1] <<- nu + w * lambda_futuro
      }
    }
    
    # Atualizar os nós dependentes e a log-probabilidade
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = TRUE)
  },
  
  methods = list(
    reset = function() {}
  )
)
