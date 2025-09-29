dynamic_sampler <- nimbleFunction(
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    dims_Y <- dim(model$Y)   # Y é [n_regions, n_times]
    n_regions <- dims_Y[1]
    n_times   <- dims_Y[2]
    w = control$w
    dims_x <- dim(model$x)   # x é [n_regions, n_times, p]
    p <- dims_x[3]
    calcNodes <- model$getDependencies(target,self = FALSE)
    targetNodes <- model$expandNodeNames(target)
    setupOutputs(
      p=p,
      n_regions=n_regions,
      n_times=n_times
    )
  },
  
  run = function() {
    declare(i, integer())
    declare(t, integer())
    declare(tt, integer())
    declare(at_buf, double())
    declare(bt_buf,double())
    declare(lambda_futuro,double())
    declare(nu,double())
    for(i in 1:n_regions) {
      shape_tmp <- model$att[i, n_times]
      rate_tmp  <- model$btt[i,n_times]
      model$lambda[i, n_times] <<- rgamma(1, shape_tmp, rate_tmp)
      
      ## --- backward recursion ---
      for(t in n_times:2) {
        at_buf <- model$at[i,t-1]
        bt_buf <- model$bt[i,t-1]
        shape_tmp <- (1 - w) * at_buf
        rate_tmp  <- bt_buf
        lambda_futuro <- model$lambda[i, t]
        nu <- rgamma(1, shape_tmp, rate_tmp)
        model$lambda[i, t-1] <<- nu + w * lambda_futuro
      }
    }
    
    ## atualizar likelihood
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = FALSE)
  },
  
  methods = list(
    reset = function() {}
  )
)
