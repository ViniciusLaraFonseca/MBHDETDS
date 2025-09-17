dynamic_sampler <- nimbleFunction(
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    ## constants passed in via control or from model
    n_regions <- control$n_regions
    n_times   <- control$n_times
    w <- control$w
    ## dependencies of the target
    calcNodes <- model$getDependencies("lambda")
    a_prev_nodes <- model$getDependencies("a_prev")
    b_prev_nodes <- model$getDependencies("b_prev")
  }
  ,
  
  run = function() {
    model$calculate(a_prev_nodes)  # Recalcula tudo que afeta a_prev
    model$calculate(b_prev_nodes)
    ## --------- SCALAR DECLARATIONS ------------
    declare(i, integer(0))
    declare(tt, integer(0))
    declare(a_local, double(0))
    declare(b_local, double(0))
    declare(lambda_smooth, double(0))
    declare(lambda_next, double(0))
    ## ------------------------------------------
    
    ## Backward smoothing
    for(i in 1:n_regions) {
      ## initialize last time point
      a_local <- model$a_prev[i, n_times]
      b_local <- model$b_prev[i, n_times]
      
      model$lambda[i, n_times] <<- rgamma(1,
                                          shape = a_local,
                                          rate  = b_local)   ## careful: RATE not SCALE
    }
    
    for(tt in (n_times-1):1) {
      for(i in 1:n_regions) {
        a_local <- model$a_prev[i, tt]
        b_local <- model$b_prev[i, tt]
        
        ## compute shape for the "innovation" term
        lambda_next <- model$lambda[i, tt+1]
        tmp_shape <- (1.0 - w) * a_local+ w * lambda_next
        tmp_rate  <- b_local
        ## ensure shape > 0 numerically (small floor)
        if(tmp_shape <= 0) tmp_shape <- 1e-10
        ## FFBS smoothing step
        lambda_smooth <- rgamma(1,
                                shape = tmp_shape,
                                rate  = tmp_rate)
        
        
        model$lambda[i, tt] <<- lambda_smooth
      }
    }
    
    ## update logProb
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  
  methods = list(reset = function() {})
)
