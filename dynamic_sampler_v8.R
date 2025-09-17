dynamic_sampler <- nimbleFunction(
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    ## constants passed in via control or from model
    beta     <- model$beta
    n_regions<- control$n_regions
    n_times   <- control$n_times
    dbeta <- as.integer(control$dbeta[1])
    Xt = control$x
    E <-  control$E
    w <- control$w
    a0 <- control$a0
    b0 <- control$b0
    count <- control$count
    ## dependencies of the target
    calcNodes <- model$getDependencies("lambda")
    at_nodes <- model$getDependencies("at")
    bt_nodes <- model$getDependencies("bt")
    att_nodes <- model$getDependencies("att")
    btt_nodes <- model$getDependencies("btt")
  }
  ,
  
  run = function() {
    model$calculate("at")  # Recalcula tudo que afeta a_prev
    model$calculate("bt")
    model$calculate("att")  # Recalcula tudo que afeta a_prev
    model$calculate("beta")
    model$calculate("epsilon")
    ## --------- SCALAR DECLARATIONS ------------
    declare(i, integer(0))
    declare(t,  integer(0))
    declare(tt, integer(0))
    
    ## ------------------------------------------
    
    ## Backward smoothing
    mab   <- array(0,c(2,n_times+1,n_regions))
    att   <- array(0,c((n_times),n_regions))
    btt   <- array(0,c((n_times),n_regions))
    at    <- array(0,c((n_times+1),n_regions))
    bt    <- array(0,c((n_times+1),n_regions))
    
    #Pred:

    
    for(i in 1:n_regions) {
      at[i,1] <- a0
      bt[i,1] <- b0
      ## initialize last time point
      for(t in 2:(n_times+1)){
        att[i,t-1] <- w*at[i,t-1]
        btt[i,t-1] <- (w*bt[i,t-1]*(exp(-(sum(Xt[i,t-1,1:dbeta]*model$beta[1:dbeta]))))/(model$epsilon[i]*E[i,t-1]))
        at[i,t] <- att[i,t-1]+(count[i,t-1])
        bt[i,t] <-w*bt[i,t-1]+(1)*exp((sum(Xt[i,t-1,1:dbeta]*model$beta[1:dbeta])))*model$epsilon[i]*E[i,t-1]
        
      }
      model$lambda[i,n_times] <<- rgamma(1,att[i,n_times],btt[i,n_times])
      for(tt in n_times:2) {
        model$lambda[i,tt-1] <<- rgamma(1,(1-w)*at[i,tt-1],bt[i,tt-1])+w*model$lambda[i,tt]
        
      }
      
    
    ## update logProb
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  }},
  
  methods = list(reset = function() {})
)
