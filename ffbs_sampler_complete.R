# Baseado em viniciuslarafonseca/mbhdetds/MBHDETDS-38c08b6f2adfe60e9a6641ba7b0fd95de6559cff/dynamic_sampler.R
ffbs_sampler_complete <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # target é o nome da variável latente no nimbleCode, ex: "lambda"
    targetAsChar <- model$expandNodeNames(target) # Ex: "lambda[1:75, 1:23]"
    targetVar <- model$getVarNames(nodes = targetAsChar[1]) # Ex: "lambda"
    
    # Verifica se o alvo é uma matriz
    targetDim <- model$getDimension(targetVar)
    if(length(targetDim) != 2) stop("FFBS sampler target must be a 2D matrix (regions x time)")
    n_regions <- targetDim[1]
    n_times   <- targetDim[2]
    
    # Obter parâmetros do control list (ESSENCIAIS)
    if(is.null(control$p)) stop("Control list must provide 'p' (number of beta covariates)")
    if(is.null(control$w)) stop("Control list must provide 'w' (discount factor)")
    if(is.null(control$a0)) stop("Control list must provide 'a0' (initial shape/a)")
    if(is.null(control$b0)) stop("Control list must provide 'b0' (initial rate/b)")
    p         <- control$p
    w         <- control$w
    a0        <- control$a0
    b0        <- control$b0
    
    # Buffers internos para os parâmetros do filtro forward
    at_buf  <- nimMatrix(nrow = n_regions, ncol = n_times + 1, init = 0, type = 'double')
    bt_buf  <- nimMatrix(nrow = n_regions, ncol = n_times + 1, init = 0, type = 'double')
    
    # Nós que dependem do alvo (lambda) e precisam ser recalculados
    # Exclui o próprio alvo (lambda)
    calcNodes <- model$getDependencies(target, self = FALSE)
    
    # Nós alvo para cópia final (nome expandido, ex: "lambda[1:75, 1:23]")
    targetNodesForCopy <- targetAsChar
    
    # setupOutputs torna essas variáveis acessíveis dentro de run()
    # Importante passar model e mvSaved
    setupOutputs(model = model,
                 mvSaved = mvSaved,
                 n_regions = n_regions, n_times = n_times, p = p,
                 w = w, a0 = a0, b0 = b0,
                 at_buf = at_buf, bt_buf = bt_buf,
                 calcNodes = calcNodes,
                 targetNodesForCopy = targetNodesForCopy,
                 targetVar = targetVar # Nome da variável alvo ("lambda")
    )
  },
  
  run = function() {
    # --- Declarações explícitas ---
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
    declare(epsilon_i, double())
    declare(E_it, double())
    declare(y_it, integer()) # Y deve ser inteiro
    
    # --- Forward Filtering (Calcula e guarda at_buf, bt_buf) ---
    for(i in 1:n_regions) {
      at_buf[i, 1] <<- a0
      bt_buf[i, 1] <<- b0
      
      # Recalcula epsilon[i] - lê gamma atual do modelo
      epsilon_i <- 1.0 - sum(model$h[i, 1:K] * model$gamma[1:K])
      epsilon_i <- max(epsilon_i, 1e-10) # Segurança
      
      for(t in 1:n_times) {
        # Passos Preditivos
        att_t <- w * at_buf[i, t]
        btt_t <- w * bt_buf[i, t]
        
        # Atualização com Observação Y[i, t] (lida do modelo)
        y_it <- model$Y[i, t]
        # Tratamento Básico de NA: Se Y for NA, usar apenas o passo preditivo
        # (Pode precisar de lógica mais sofisticada dependendo do modelo)
        if(is.na(y_it)) {
          at_buf[i, t+1] <<- att_t
        } else {
          at_buf[i, t+1] <<- att_t + y_it
        }
        
        # Calcula inprod(beta, x) lendo beta e x do modelo
        prod_val <- 0.0
        for(k in 1:p) { prod_val <- prod_val + model$x[i, t, k] * model$beta[k] }
        
        E_it <- model$E[i, t] # Lê E do modelo
        # Atualiza bt (taxa/rate)
        bt_buf[i, t+1] <<- btt_t + E_it * epsilon_i * exp(prod_val)
        
        # Segurança numérica para bt_buf
        if(bt_buf[i, t+1] <= 1e-10) bt_buf[i, t+1] <<- 1e-10
      }
    }
    
    # --- Backward Sampling (Amostra lambda[1:T] e atualiza o modelo) ---
    # Matriz temporária para guardar a nova trajetória
    lambda_new <- nimMatrix(nrow=n_regions, ncol=n_times, init=FALSE)
    
    for(i in 1:n_regions) {
      # Amostra o último ponto (T) usando at_buf/bt_buf[T+1]
      shape_T <- max(at_buf[i, n_times + 1], 1e-10)
      rate_T  <- max(bt_buf[i, n_times + 1], 1e-10)
      lambda_new[i, n_times] <- rgamma(1, shape = shape_T, rate = rate_T)
      
      # Loop para trás (T-1 a 1)
      if (n_times > 1) {
        for(tt in n_times:2) { # tt vai de T até 2
          lambda_futuro <- lambda_new[i, tt] # Usa o valor recém-amostrado
          
          # Parâmetros para nu (usa at_buf/bt_buf no índice tt, que corresponde ao fim do tempo tt-1)
          shape_nu <- max((1 - w) * at_buf[i, tt], 1e-10)
          rate_nu  <- max(bt_buf[i, tt], 1e-10)
          
          nu <- rgamma(1, shape = shape_nu, rate = rate_nu)
          
          # Reconstrói lambda no tempo tt-1
          lambda_new[i, tt-1] <- nu + w * lambda_futuro
        }
      }
    }
    
    # --- Atualizar Modelo e Copiar Estado ---
    # Atribui a matriz completa da nova trajetória ao nó alvo no modelo
    values(model, targetVar) <<- lambda_new
    
    # Recalcula nós dependentes (Y, mu, theta, logProbs...) no modelo
    # Isso atualiza a verossimilhança com base na nova trajetória lambda
    lp <- model$calculate(calcNodes) # Guarda o logProb calculado
    
    # Copia o estado atualizado (incluindo lambda E dependências) para mvSaved
    # Copia calcNodes primeiro (com logProb)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    # Copia o nó alvo amostrado (lambda) depois (sem logProb, pois já está implícito no lp total)
    copy(from = model, to = mvSaved, row = 1, nodes = targetNodesForCopy, logProb = FALSE)
    
    # Opcional: Retornar o logProb calculado pode ser útil para alguns diagnósticos
    # return(lp) # Descomente se precisar
    
  },
  methods = list( reset = function () {} )
)
print("Amostrador FFBS completo (ffbs_sampler_complete) definido.")