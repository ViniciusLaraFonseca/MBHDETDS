# -------------------------------------------------------------------
# AMOSTRADOR FFBS CUSTOMIZADO (NIMBLEFUNCTION)
# PARA O MODELO lambda_it (POR ÁREA-TEMPO)
# -------------------------------------------------------------------

ffbs_sampler_Lambda_it <- nimbleFunction(
  
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    
    # --- Extrair dimensões e parâmetros ---
    n_regions <- control$n_regions
    n_times   <- control$n_times
    p         <- control$p
    w         <- control$w
    a0        <- control$a0
    b0        <- control$b0
    
    # --- Buffers 2D (Matrizes) ---
    at_buf  <- nimNumeric(nrow = n_regions, ncol = n_times + 1, value = 0)
    bt_buf  <- nimNumeric(nrow = n_regions, ncol = n_times + 1, value = 0)
    
    # --- Nós dependentes ---
    calcNodes   <- model$getDependencies(target, self = FALSE)
    targetNodes <- model$expandNodeNames(target) # target é 'lambda'
    
    # --- Retornar buffers ---
    setupOutputs(at_buf, bt_buf) 
  },
  
  run = function() {
    
    # --- Declarações ---
    declare(i, integer())
    declare(t, integer())
    declare(k, integer())
    declare(prod_val, double())
    declare(g_it, double())
    declare(att_t, double())
    declare(btt_t, double())
    declare(shape_tmp, double())
    declare(rate_tmp, double())
    declare(lambda_futuro, double())
    declare(nu, double())
    
    # --- Loop principal sobre as ÁREAS ---
    for(i in 1:n_regions) {
      
      # ================================================================
      # PASSO 1 — FORWARD FILTERING (para a área i)
      # ================================================================
      
      at_buf[i, 1] <<- a0
      bt_buf[i, 1] <<- b0
      
      for(t in 1:n_times) {
        
        # 1.1 — Preditiva (para i, t)
        att_t <- w * at_buf[i, t]
        btt_t <- w * bt_buf[i, t]
        
        # 1.2 — Componente espacial g_it
        prod_val <- 0.0
        for(k in 1:p) {
          prod_val <- prod_val + model$x[i, t, k] * model$beta[k]
        }
        g_it <- model$E[i, t] * model$epsilon[i] * exp(prod_val)
        
        # 1.3 — Atualização (para i, t)
        at_buf[i, t + 1] <<- att_t + model$Y[i, t] # Y[i,t]
        bt_buf[i, t + 1] <<- btt_t + g_it          # g_it
      }
      
      # ================================================================
      # PASSO 2 — BACKWARD SAMPLING (para a área i)
      # ================================================================
      
      # 2.1 — Amostrar λ_{i,T}
      shape_tmp <- at_buf[i, n_times + 1]
      rate_tmp  <- bt_buf[i, n_times + 1]
      model$lambda[i, n_times] <<- rgamma(1, shape = shape_tmp, rate = rate_tmp)
      
      # 2.2 — Amostrar λ_{i,t}, t = T-1 ... 1
      for(t_idx in (n_times - 1):1) {
        lambda_futuro <- model$lambda[i, t_idx + 1] # λ_{i, t+1}
        
        shape_tmp <- (1 - w) * at_buf[i, t_idx + 1] # (1-w) * a_it
        rate_tmp  <- bt_buf[i, t_idx + 1]           # b_it
        
        nu <- rgamma(1, shape = shape_tmp, rate = rate_tmp)
        
        model$lambda[i, t_idx] <<- nu + w * lambda_futuro
      }
      
    } # --- Fim do loop sobre as ÁREAS ---
    
    # ================================================================
    # PASSO 3 — ATUALIZAÇÃO (depois que TODOS os lambdas foram amostrados)
    # ================================================================
    
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, 
         nodes = targetNodes, logProb = TRUE)
  },
  
  methods = list(
    reset = function() {}
  )
)