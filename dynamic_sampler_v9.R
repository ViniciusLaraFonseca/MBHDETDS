
if(exists("nimbleFunction", mode = "function")) rm(nimbleFunction)
if(exists("dynamic_sampler")) rm(dynamic_sampler)
if(!require(nimble)){ install.packages("nimble", repos = "http://cran.r-project.org"); require(nimble)}




dynamic_sampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    n_regions <- as.integer(control$n_regions)[1]
    n_times   <- as.integer(control$n_times)[1]
    dbeta     <- as.integer(control$dbeta)[1]
    w         <- control$w
    a0        <- control$a0
    b0        <- control$b0
    
    ## buffers 2D via nimMatrix
    att_buf <- nimMatrix(init = 0, nrow = n_regions, ncol = n_times)
    btt_buf <- nimMatrix(init = 0, nrow = n_regions, ncol = n_times)
    at_buf  <- nimMatrix(init = 0, nrow = n_regions, ncol = n_times + 1)
    bt_buf  <- nimMatrix(init = 0, nrow = n_regions, ncol = n_times + 1)
  },
  
  run = function() {
    declare(i, integer(0))
    declare(t, integer(0))
    declare(tt, integer(0))
    declare(k, integer(0))
    
    # Ensure these are accessible
    declare(n_regions, integer(0))
    declare(n_times, integer(0))
    declare(dbeta, integer(0))
    declare(w, double(0))
    declare(a0, double(0))
    declare(b0, double(0))
    
    model$calculate("epsilon")
    calcNodes <- model$getDependencies("lambda")
    model$calculate(calcNodes)
    
    for(i in 1:n_regions) {
      at_buf[i, 1] <<- a0
      bt_buf[i, 1] <<- b0
      
      prod_val <- 0
      for(k in 1:dbeta) {
        prod_val <- prod_val + model$x[i, 1, k] * model$beta[k]
      }
      
      for(t in 2:(n_times + 1)) {
        att_buf[i, t-1] <<- w * at_buf[i, t-1]
        btt_buf[i, t-1] <<- w * bt_buf[i, t-1]
        at_buf[i, t]    <<- att_buf[i, t-1] + model$count[i, t-1]
        
        prod_val <- 0
        for(k in 1:dbeta) {
          prod_val <- prod_val + model$x[i, t-1, k] * model$beta[k]
        }
        bt_buf[i, t] <<- btt_buf[i, t-1] +
          model$epsilon[i] * model$E[i, t-1] * exp(prod_val)
      }
      
      model$lambda[i, n_times] <<-
        rgamma(1, at_buf[i, n_times], btt_buf[i, n_times])
      
      for(tt in n_times:2) {
        model$lambda[i, tt-1] <<-
          rgamma(1,
                 (1 - w) * at_buf[i, tt-1],
                 bt_buf[i, tt-1]
          ) + w * model$lambda[i, tt]
      }
      
      copy(from = model, to = mvSaved, row = 1,
           nodes = calcNodes, logProb = TRUE)
    }
  },
  
  methods = list(
    reset = function() {}
  )
)

registerSampler(name = 'dynamic', sampler = dynamic_sampler)

