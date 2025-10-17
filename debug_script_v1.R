# --- PASSO 0: CARREGAR PACOTES E GARANTIR REPRODUTIBILIDADE ---
if (!require(nimble)) {
  install.packages("nimble")
}
if (!require(coda)) {
  install.packages("coda")
}
library(nimble)
library(coda)

set.seed(123)

cat("--- Início do Script de Depuração (Para Gerar Prints) ---\n\n")

# --- PASSO 1: DADOS REDUZIDOS E PARÂMETROS ---
cat("--- PASSO 1: Configurando ambiente de teste reduzido ---\n")

n_regions_debug <- 1; n_times_debug <- 5; p_debug <- 3; K_debug <- 4
beta_initial_values <- c(-0.1, 0.2, -0.3)
gamma_fixed <- c(0.05, 0.10, 0.10, 0.15)
w_debug <- 0.9; a0_debug <- 1.0; b0_debug <- 1.0

# Gerar dados fictícios, incluindo a normalização de E
E_raw <- matrix(runif(n_regions_debug * n_times_debug, 100, 200), nrow = n_regions_debug)
mean_E <- mean(E_raw)
E_debug <- E_raw / mean_E # Normalizar dividindo pela média
cat("Matriz de Offset 'E' normalizada (média 1):\n"); print(E_debug); cat("\n")

x_debug <- array(rnorm(n_regions_debug * n_times_debug * p_debug), dim = c(n_regions_debug, n_times_debug, p_debug))
hAI_debug <- matrix(c(1, 0, 0, 0), nrow = n_regions_debug)
epsilon_fixed <- 1 - as.numeric(hAI_debug %*% gamma_fixed)
lambda_true_sim <- rgamma(n_times_debug, shape = 2, rate = 0.5)
mu_true_sim <- E_debug[1,] * epsilon_fixed * exp(x_debug[1,,] %*% beta_initial_values) * lambda_true_sim
y_ini_debug <- matrix(rpois(n_times_debug, lambda = mu_true_sim), nrow = n_regions_debug)

cat("Dados de Y (observações) para a depuração:\n"); print(y_ini_debug); cat("\n")

# --- PASSO 2: DYNAMIC_SAMPLER INSTRUMENTADO ---
cat("--- PASSO 2: Definindo o amostrador 'dynamic_sampler_debug' com impressões detalhadas ---\n")

dynamic_sampler_debug <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    dims_Y <- dim(model$Y); n_regions <- dims_Y[1]; n_times <- dims_Y[2]
    p <- dim(model$x)[3]; w <- control$w; a0 <- control$a0; b0 <- control$b0
    at_buf  <- nimMatrix(nrow = n_regions, ncol = n_times + 1, init = 0, type = 'double')
    bt_buf  <- nimMatrix(nrow = n_regions, ncol = n_times + 1, init = 0, type = 'double')
    calcNodes   <- model$getDependencies(target, self = FALSE)
    targetNodes <- model$expandNodeNames(target)
    setupOutputs(n_regions, n_times, p, w, a0, b0, at_buf, bt_buf, calcNodes, targetNodes)
  },
  run = function() {
    declare(i, integer()); declare(t, integer()); declare(tt, integer()); declare(k, integer())
    declare(prod_val, double()); declare(att_t, double()); declare(btt_t, double())
    declare(shape_tmp, double()); declare(rate_tmp, double()); declare(lambda_futuro, double())
    declare(nu, double())
    
    for(i in 1:n_regions) {
      at_buf[i, 1] <<- a0
      bt_buf[i, 1] <<- b0
      
      print(paste0("--- INICIANDO FORWARD FILTERING (ITERAÇÃO MCMC) PARA REGIÃO ", i, " ---"))
      
      for(t in 1:n_times) {
        att_t <- w * at_buf[i, t]
        btt_t <- w * bt_buf[i, t]
        
        at_buf[i, t+1] <<- att_t + model$Y[i, t]
        
        prod_val <- 0
        for(k in 1:p) {
          prod_val <- prod_val + model$x[i, t, k] * model$beta[k]
        }
        
        bt_buf[i, t+1] <<- btt_t + model$E[i, t] * model$epsilon[i] * exp(prod_val)
        
        print(paste0("Forward, t=", t, ": at_post=", round(at_buf[i, t+1], 4), ", bt_post=", round(bt_buf[i, t+1], 4)))
      }
    }
    
    for(i in 1:n_regions) {
      print(paste0("--- INICIANDO BACKWARD SAMPLING (ITERAÇÃO MCMC) PARA REGIÃO ", i, " ---"))
      
      shape_tmp_final <- at_buf[i, n_times + 1]
      rate_tmp_final  <- bt_buf[i, n_times + 1]
      model$lambda[i, n_times] <<- rgamma(1, shape = shape_tmp_final, rate = rate_tmp_final)
      print(paste0("Backward, t=", n_times, ": lambda=", round(model$lambda[i, n_times], 6), 
                   " (de Gamma(shape=", round(shape_tmp_final, 4), ", rate=", round(rate_tmp_final, 4), "))"))
      
      for(tt in n_times:2) { 
        lambda_futuro <- model$lambda[i, tt]
        
        shape_tmp <- (1 - w) * at_buf[i, tt] 
        rate_tmp  <- bt_buf[i, tt]
        
        if(shape_tmp <= 0 | rate_tmp <= 0) {
          print(paste0("!!! ERRO: Parâmetros da Gamma inválidos em t=", tt-1, "!!!"))
        }
        
        nu <- rgamma(1, shape = shape_tmp, rate = rate_tmp)
        model$lambda[i, tt-1] <<- nu + w * lambda_futuro
        
        print(paste0("Backward, t=", tt-1, ": nu=", format(nu, scientific = TRUE, digits = 4), 
                     " (de Gamma(shape=", round(shape_tmp, 4), ", rate=", round(rate_tmp, 4), "))",
                     " | lambda=", round(model$lambda[i, tt-1], 6)))
      }
    }
    
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = TRUE)
  },
  methods = list(reset = function() {})
)

# --- PASSO 3: CONFIGURAR E EXECUTAR O MODELO DE DEPURAÇÃO (NÃO COMPILADO) ---
cat("\n--- PASSO 3: Configurando e executando o modelo NIMBLE em R ---\n")
code_debug <- nimbleCode({
  for (j in 1:p) { beta[j] ~ dnorm(0, sd = 1.5) } 
  for (j in 1:K) { gamma[j] ~ dunif(0, 1) } 
  for (i in 1:n_regions) { epsilon[i] <- 1 - sum(h[i, 1:K] * gamma[1:K]) }
  for (i in 1:n_regions) {
    for(t in 1:n_times){  
      lambda[i, t] ~ dgamma(1, 1)
      log_theta[i, t] <- log(lambda[i, t]) + inprod(beta[1:p], x[i, t, 1:p])
      theta[i, t] <- exp(log_theta[i, t])
      mu[i,t] <- E[i, t] * epsilon[i] * theta[i, t]
      Y[i, t] ~ dpois(mu[i,t])
    }
  }
})
constants_debug <- list(n_regions = n_regions_debug, n_times = n_times_debug, p = p_debug, K = K_debug, h = hAI_debug)
data_debug <- list(Y = y_ini_debug, E = E_debug, x = x_debug, gamma = gamma_fixed) 
inits_debug <- list(lambda = matrix(1, nrow = n_regions_debug, ncol = n_times_debug), beta = beta_initial_values)

model_debug <- nimbleModel(code_debug, constants = constants_debug, data = data_debug, inits = inits_debug, check = FALSE)

conf_debug <- configureMCMC(model_debug, nodes = NULL)
conf_debug$addSampler(target = "lambda", type = dynamic_sampler_debug, control = list(w = w_debug, a0 = a0_debug, b0 = b0_debug))
conf_debug$addSampler(target = "beta", type = "RW_block", control = list(adaptInterval = 200))
conf_debug$monitors <- c("lambda", "beta")

# Construir a versão R do MCMC (não compilada)
Rmcmc_debug <- buildMCMC(conf_debug)

cat("\n--- EXECUTANDO MCMC (VERSÃO R) POR 5 ITERAÇÕES PARA GERAR PRINTS ---\n\n")
# Executar por poucas iterações para analisar a saída detalhada
# Nota: NÃO compilamos o MCMC para que os prints funcionem na consola
runMCMC(Rmcmc_debug, niter = 5, nburnin = 0, nchains = 1, summary = FALSE)
cat("\n--- EXECUÇÃO DOS PRINTS CONCLUÍDA ---\n")