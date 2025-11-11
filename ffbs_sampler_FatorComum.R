# -------------------------------------------------------------------
# AMOSTRADOR FFBS CUSTOMIZADO (NIMBLEFUNCTION)
# PARA O MODELO DE FATOR DINÂMICO COMUM (LAMBDA_T UNIVARIADO)
# -------------------------------------------------------------------
# OBS: NÃO carregamos library(nimble) aqui — o script principal deve
# garantir que 'nimble' esteja anexado (library(nimble)) antes do source().
# -------------------------------------------------------------------

ffbs_sampler_FatorComum <- nimbleFunction(
  
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    
    n_regions <- control$n_regions
    n_times   <- control$n_times
    p         <- control$p
    w         <- control$w
    a0        <- control$a0
    b0        <- control$b0
    
    # Buffers corretos (vetores R normais)
    at_buf <- numeric(n_times + 1)
    bt_buf <- numeric(n_times + 1)
    
    calcNodes   <- model$getDependencies(target, self = FALSE)
    targetNodes <- model$expandNodeNames(target)
  },
  
  run = function() {
    
    # --- Declaração explícita de tipos ---
    declare(i, integer())
    declare(t, integer())
    declare(k, integer())
    declare(prod_val, double())
    declare(g_it, double())
    declare(sum_Y_t, double())
    declare(sum_g_t, double())
    declare(att_t, double())
    declare(btt_t, double())
    declare(shape_tmp, double())
    declare(rate_tmp, double())
    declare(lambda_futuro, double())
    declare(nu, double())
    
    # ================================================================
    # PASSO 1 — FORWARD FILTERING
    # ================================================================
    
    # Inicialização no tempo 0
    at_buf[1] <<- a0
    bt_buf[1] <<- b0
    
    for(t in 1:n_times) {
      
      # 1.1 — Passo Preditivo
      att_t <- w * at_buf[t]
      btt_t <- w * bt_buf[t]
      
      # 1.2 — Agregação
      sum_Y_t <- 0.0
      sum_g_t <- 0.0
      
      for(i in 1:n_regions) {
        prod_val <- 0.0
        for(k in 1:p) {
          prod_val <- prod_val + model$x[i, t, k] * model$beta[k]
        }
        
        g_it <- model$E[i, t] * model$epsilon[i] * exp(prod_val)
        
        sum_g_t <- sum_g_t + g_it
        sum_Y_t <- sum_Y_t + model$Y[i, t]
      }
      
      # 1.3 — Atualização
      at_buf[t + 1] <<- att_t + sum_Y_t
      bt_buf[t + 1] <<- btt_t + sum_g_t
    }
    
    # ================================================================
    # PASSO 2 — BACKWARD SAMPLING
    # ================================================================
    
    # 2.1 — λ_T
    shape_tmp <- at_buf[n_times + 1]
    rate_tmp  <- bt_buf[n_times + 1]
    model$lambda[n_times] <<- rgamma(1, shape = shape_tmp, rate = rate_tmp)
    
    # 2.2 — λ_{t}, t = T-1 ... 1
    for(t_idx in (n_times - 1):1) {
      
      lambda_futuro <- model$lambda[t_idx + 1]
      
      shape_tmp <- (1 - w) * at_buf[t_idx + 1]
      rate_tmp  <- bt_buf[t_idx + 1]
      
      nu <- rgamma(1, shape = shape_tmp, rate = rate_tmp)
      
      model$lambda[t_idx] <<- nu + w * lambda_futuro
    }
    
    # ================================================================
    # PASSO 3 — ATUALIZAÇÃO DO LOG-PROB
    # ================================================================
    
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, 
         nodes = targetNodes, logProb = TRUE)
  },
  
  methods = list(
    reset = function() {}
  )
)
