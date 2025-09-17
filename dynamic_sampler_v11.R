dynamic_sampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## checagens de controle
    if(is.null(control$n_regions)) stop("n_regions não definido no controle")
    if(is.null(control$n_times))   stop("n_times não definido no controle")
    if(is.null(control$dbeta))     stop("dbeta não definido no controle")
    if(is.null(control$w))         stop("w não definido no controle")
    if(is.null(control$a0))        stop("a0 não definido no controle")
    if(is.null(control$b0))        stop("b0 não definido no controle")
    n_regions <- as.integer(control$n_regions)[1]
    n_times   <- as.integer(control$n_times)[1]
    dbeta     <- as.integer(control$dbeta)[1]
    w         <- control$w
    a0        <- constants$a0
    b0        <- constants$b0
    ## buffers 2D via nimMatrix
    att_buf <- nimMatrix(init = 0, nrow = n_regions, ncol = n_times)
    btt_buf <- nimMatrix(init = 0, nrow = n_regions, ncol = n_times)
    at_buf  <- nimMatrix(init = 0, nrow = n_regions, ncol = n_times + 1)
    bt_buf  <- nimMatrix(init = 0, nrow = n_regions, ncol = n_times + 1)
    ## nós dependentes do target
    calcNodes <- model$getDependencies(target)
    ## tudo que o run() precisa tem que ser devolvido aqui
    setupOutputs(
      
        n_regions = n_regions,
        n_times   = n_times,
        dbeta     = dbeta,
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
    declare(i,  integer(0))
    declare(t,  integer(0))
    declare(tt, integer(0))
    declare(k,  integer(0))
    model$calculate(nodes = "epsilon")
    for(i in 1:n_regions) {
      at_buf[i, 1] <<- a0
      bt_buf[i, 1] <<- b0
      prod_val <- 0
      for(k in 1:dbeta) {
        prod_val <- prod_val + model$x[i, 1, k] * model$beta[k]
      }
      for(t in 2:(n_times + 1)) {
        att_buf[i, t-1] <<- w * at_buf[i, t-1]
        btt_buf[i, t-1] <<- w * bt_buf[i, t-1]/(model$epsilon[i] * model$E[i, t-1] * exp(prod_val))
        at_buf[i, t]    <<- att_buf[i, t-1] + model$count[i, t-1]
        prod_val <- 0
        for(k in 1:dbeta) {
          prod_val <- prod_val + model$x[i, t-1, k] * model$beta[k]
        }
        bt_buf[i, t] <<- w*bt_buf[i, t-1] +
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
    }
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1,
         nodes = calcNodes, logProb = TRUE)
  },
  methods = list(
    reset = function() {}
  )
)
