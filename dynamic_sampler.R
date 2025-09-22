dynamic_sampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## checagens de controle
    n_regions <- as.integer(constants$n_regions)[1]
    n_times   <- as.integer(constants$n_times)[1]
    p         <- as.integer(constants$p)[1]
    w         <- constants$w
    a0        <- constants$a0
    b0        <- constants$b0
    
    ## buffers 2D via nimMatrix com tipo explícito
    att_buf <- nimMatrix(nrow = n_regions, ncol = n_times, init = 0, type = 'double')
    btt_buf <- nimMatrix(nrow = n_regions, ncol = n_times, init = 0, type = 'double')
    at_buf  <- nimMatrix(nrow = n_regions, ncol = (n_times+1), init = 0, type = 'double')
    bt_buf  <- nimMatrix(nrow = n_regions, ncol = (n_times+1), init = 0, type = 'double')
    
    ## nós dependentes do target
    calcNodes <- control$calcNodes
    
    ## declarar variáveis explicitamente
    
    
    ## devolver objetos para run()
    setupOutputs(
      n_regions = n_regions,
      n_times   = n_times,
      p         = p,
      w         = w,
      a0        = a0,
      b0        = b0,
      att_buf   = att_buf,
      btt_buf   = btt_buf,
      at_buf    = at_buf,
      bt_buf    = bt_buf,
      calcNodes = calcNodes
    )
  },
  
  run = function() {
    declare(i, integer())
    declare(t, integer())
    declare(tt, integer())
    declare(k, integer())
    declare(prod_val, double())
    
    for(i in 1:n_regions) {
      at_buf[i, 1] <<- a0
      bt_buf[i, 1] <<- b0
      
      ## Inicializar prod_val
      prod_val <- 0
      for(k in 1:p) {
        prod_val <- prod_val + model$x[i, 1, k] * model$beta[k]
      }
      
      for(t in 2:(n_times + 1)) {
        att_buf[i, t-1] <<- w * at_buf[i, t-1]
        btt_buf[i, t-1] <<- w * bt_buf[i, t-1] / (model$epsilon[i] * model$E[i, t-1] * exp(prod_val))
        at_buf[i, t]    <<- att_buf[i, t-1] + model$count[i, t-1]
        
        ## recalcular prod_val
        prod_val <- 0
        for(k in 1:p) {
          prod_val <- prod_val + model$x[i, t-1, k] * model$beta[k]
        }
        
        bt_buf[i, t] <<- w * bt_buf[i, t-1] +
          model$epsilon[i] * model$E[i, t-1] * exp(prod_val)
      }
      
      ## amostragem lambda (ultimo tempo)
      model$lambda[i, n_times] <<- rgamma(1, at_buf[i, n_times], btt_buf[i, n_times])
      
      ## backward recursion
      for(tt in seq(n_times, 2, by = -1)) {
        model$lambda[i, tt-1] <<- 
          rgamma(1, (1 - w) * at_buf[i, tt-1], bt_buf[i, tt-1]) +
          w * model$lambda[i, tt]
      }
    }
    
    ## calcular logProb
    model$calculate(calcNodes)
    
    ## copiar para mvSaved
    copy(from = model, to = mvSaved, row = 1,
         nodes = calcNodes, logProb = TRUE)
  },
  
  methods = list(
    reset = function() {}
  )
)
