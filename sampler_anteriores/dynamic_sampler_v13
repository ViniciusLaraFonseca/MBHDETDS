dynamic_sampler <- nimbleFunction(
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    dims_Y <- dim(model$Y)   # Y é [n_regions, n_times]
    n_regions <- dims_Y[1]
    n_times   <- dims_Y[2]
    
    dims_x <- dim(model$x)   # x é [n_regions, n_times, p]
    p <- dims_x[3]
    
    # --- hiperparâmetros vindos do control ---
    w  <- control$w
    a0 <- control$a0
    b0 <- control$b0
    
    # --- buffers para recursão ---
    att_buf <- nimMatrix(nrow = n_regions, ncol = n_times, init = 0, type = 'double')
    btt_buf <- nimMatrix(nrow = n_regions, ncol = n_times, init = 0, type = 'double')
    at_buf  <- nimMatrix(nrow = n_regions, ncol = n_times+1, init = 0, type = 'double')
    bt_buf  <- nimMatrix(nrow = n_regions, ncol = n_times+1, init = 0, type = 'double')
    
    # --- dependências ---
    calcNodes <<- model$getDependencies(target)
    setupOutputs(
      n_regions = n_regions,
      n_times   = n_times,
      p         = p,
      w  = w, a0 = a0, b0 = b0,
      att_buf = att_buf,
      btt_buf = btt_buf,
      at_buf  = at_buf,
      bt_buf  = bt_buf
      
    )

  },
  
  run = function() {
    declare(i, integer())
    declare(t, integer())
    declare(tt, integer())
    declare(k, integer())
    declare(prod_val, double())
    declare(tmp, double())
    declare(shape_tmp, double())
    declare(rate_tmp, double())
    declare(E_it,double())
    declare(count_it,double())
    declare(epsilon_i,double())
    
    for(i in 1:n_regions) {
      at_buf[i,1] <<- a0
      bt_buf[i,1] <<- b0
      epsilon_i <- model$epsilon[i]
      ## --- forward recursion ---
      for(t in 2:(n_times+1)) {
        prod_val <- 0
        for(k in 1:p) {
          tmp <- model$x[i, t-1, k]
          prod_val <- prod_val + tmp * model$beta[k]
        }
        count_it <- model$count[i,t-1]
        E_it <- model$E[i,t-1]
        att_buf[i, t-1] <<- w * at_buf[i, t-1]
        btt_buf[i, t-1] <<- w * bt_buf[i, t-1] / (epsilon_i * E_it * exp(prod_val))
        
        at_buf[i, t] <<- att_buf[i, t-1] + count_it
        bt_buf[i, t] <<- btt_buf[i, t-1] + 1
      }
      
      ## --- amostragem lambda no último tempo ---
      shape_tmp <- att_buf[i, n_times]
      rate_tmp  <- btt_buf[i, n_times]
      model$lambda[i, n_times] <<- rgamma(1, shape_tmp, rate_tmp)
      
      ## --- backward recursion ---
      for(tt in n_times:2) {
        shape_tmp <- (1 - w) * at_buf[i, tt-1]
        rate_tmp  <- bt_buf[i, tt-1]
        model$lambda[i, tt-1] <<- rgamma(1, shape_tmp, rate_tmp) + w * model$lambda[i, tt]
      }
    }
    
    ## atualizar likelihood
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  
  methods = list(
    reset = function() {}
  )
)

