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

cat("--- Início do Script de Depuração Final (Block Sampler AF_slice) ---\n\n")

# --- PASSO 1: DADOS REDUZIDOS E PARÂMETROS ---
n_regions_debug <- 1; n_times_debug <- 20; p_debug <- 3; K_debug <- 4
beta_initial_values <- c(0.5, 1, -0.3)
gamma_fixed <- c(0.05, 0.10, 0.10, 0.15)
w_debug <- 0.9; # Voltamos para w=0.9 para testar o modelo original
a0_debug <- 1.0; b0_debug <- 1.0
E_debug_at <- matrix(runif(n_regions_debug * (n_times_debug+1), 100, 200), nrow = n_regions_debug)
E_debug <- matrix(E_debug_at[,2:(n_times_debug+1)],nrow=n_regions_debug)
x_debug_at <- array(rnorm(n_regions_debug * (n_times_debug+1) * p_debug), dim = c(n_regions_debug, (n_times_debug+1), p_debug))
x_debug <- array(x_debug_at[,2:(n_times_debug+1),],dim = c(n_regions_debug, (n_times_debug+1), p_debug))
hAI_debug <- matrix(c(1, 0, 0, 0), nrow = n_regions_debug)
epsilon_fixed <- 1 - as.numeric(hAI_debug %*% gamma_fixed)
#EFEITO DINAMICO 
lambda0 <- rgamma(n_times_debug, shape = a0_debug, rate = b0_debug)
mu_true0<- E_debug[1,1] * epsilon_fixed * exp(c(x_debug[1,1,] %*% beta_initial_values)) * lambda0
y_ini_debug0 <- matrix(rpois(1, lambda = mu_true0), nrow = n_regions_debug)
lambda_at_fwd <- matrix(NA, ncol = n_times_debug, nrow = n_regions_debug)
y_at <- matrix(NA, nrow = n_regions_debug, ncol = n_times_debug + 1)
y_at[, 1] <- y_ini_debug0  # Inicia com a média observada por região

att_ini <- matrix(NA, nrow = n_regions_debug, ncol = n_times_debug)
btt_ini <- matrix(NA, nrow = n_regions_debug, ncol = n_times_debug)
at_ini  <- matrix(NA, nrow = n_regions_debug, ncol = n_times_debug + 1)
bt_ini  <- matrix(NA, nrow = n_regions_debug, ncol = n_times_debug + 1)

# Simulação Forward
for(i in 1:n_regions_debug) {
  at_ini[i, 1] <- a0_debug
  bt_ini[i, 1] <- b0_debug
  for(t in 2:(n_times_debug + 1)) {
    # Passos de filtragem (forward)
    att_ini[i, t-1] <- w_debug * at_ini[i, t-1]
    btt_ini[i, t-1] <- w_debug * bt_ini[i, t-1]
    at_ini[i, t]    <- att_ini[i, t-1] + y_at[i, t-1]
    
    prod_val <- sum(x_debug_at[i, t, ] * beta_initial_values)
    bt_ini[i, t]    <- btt_ini[i, t-1] + E_debug_at[i, t] * epsilon_fixed[i] * exp(prod_val)
    
    # Amostra um valor de lambda a partir da distribuição filtrada (forward)
    lambda_at_fwd[i,t-1] <- rgamma(1, shape = at_ini[i, t], rate = bt_ini[i, t])
    
    # Simula um novo valor de Y para alimentar o próximo passo do filtro
    mu_at <- lambda_at_fwd[i, t-1] * exp(prod_val) * epsilon_fixed[i] * E_debug_at[i, t]
    y_at[i, t] <- rpois(1, mu_at)
  }
}
y_ini_debug <- y_at[,2:21, drop = FALSE]

# --- PASSO 2: DYNAMIC_SAMPLER (NÃO PRECISA DE ALTERAÇÕES) ---
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
      at_buf[i, 1] <<- a0; bt_buf[i, 1] <<- b0
      for(t in 1:n_times) {
        att_t <- w * at_buf[i, t]; btt_t <- w * bt_buf[i, t]
        at_buf[i, t+1] <<- att_t + model$Y[i, t]
        prod_val <- 0
        for(k in 1:p) { prod_val <- prod_val + model$x[i, t, k] * model$beta[k] }
        bt_buf[i, t+1] <<- btt_t + model$E[i, t] * model$epsilon[i] * exp(prod_val)
      }
    }
    for(i in 1:n_regions) {
      shape_tmp_final <- at_buf[i, n_times + 1]; rate_tmp_final  <- bt_buf[i, n_times + 1]
      model$lambda[i, n_times] <<- rgamma(1, shape = shape_tmp_final, rate = rate_tmp_final)
      for(tt in n_times:2) {
        lambda_futuro <- model$lambda[i, tt]
        shape_tmp <- (1 - w) * at_buf[i, tt]; rate_tmp  <- bt_buf[i, tt]
        nu <- rgamma(1, shape = shape_tmp, rate = rate_tmp)
        model$lambda[i, tt-1] <<- nu + w * lambda_futuro
      }
    }
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = TRUE)
  },
  methods = list(reset = function() {})
)

# --- PASSO 3: CONFIGURAR E EXECUTAR O MODELO COM AF_slice SAMPLER ---
cat("\n--- PASSO 3: Configurando e compilando com o amostrador AF_slice para beta ---\n")
code_debug <- nimbleCode({
  for (j in 1:p) { beta[j] ~ dnorm(0, sd = 1.5) } 
  for (j in 1:K) { gamma[j] ~ dunif(0, 1) } 
  for (i in 1:n_regions) { epsilon[i] <- 1 - sum(h[i, 1:K] * gamma[1:K]) }
  for (i in 1:n_regions) {
    for(t in 1:n_times){  
      lambda[i, t] ~ dgamma(1, 1)
      log_theta[i, t] <- log(lambda[i, t]) + inprod(x[i, t, 1:p],beta[1:p])
      theta[i, t] <- exp(log_theta[i, t])
      mu[i,t] <- E[i, t] * epsilon[i] * theta[i, t]
      Y[i, t] ~ dpois(mu[i,t])
    }
  }
})
constants_debug <- list(n_regions = n_regions_debug, n_times = n_times_debug, p = p_debug, K = K_debug, h = hAI_debug)
data_debug <- list(Y = y_ini_debug, E = E_debug, x = x_debug)
inits_debug <- list(lambda = matrix(1, nrow = n_regions_debug, ncol = n_times_debug), beta = beta_initial_values,gamma = c(0.05, 0.10, 0.10, 0.15))

model_debug <- nimbleModel(code_debug, constants = constants_debug, data = data_debug, inits = inits_debug, check = FALSE)
Cmodel_debug <- compileNimble(model_debug)

conf_debug <- configureMCMC(Cmodel_debug, nodes = NULL)
conf_debug$addSampler(target = "lambda", type = dynamic_sampler_debug, control = list(w = w_debug, a0 = a0_debug, b0 = b0_debug))

# ALTERAÇÃO FINAL: Usar o amostrador AF_slice para o vetor beta inteiro
conf_debug$removeSamplers("beta") 
conf_debug$addSampler(target = "beta", type = "AF_slice")

conf_debug$monitors <- c("lambda", "beta")

cat("\nConfiguração do MCMC com o amostrador AF_slice:\n")
conf_debug$printSamplers()

Rmcmc_debug <- buildMCMC(conf_debug)
Cmcmc_debug <- compileNimble(Rmcmc_debug, project = Cmodel_debug, resetFunctions = TRUE)
cat("\n--- Compilação concluída com sucesso! ---\n")

cat("\n--- Executando MCMC (VERSÃO C++) por 10000 iterações ---\n")
mcmc.out_compiled <- runMCMC(Cmcmc_debug, niter = 10000, nburnin = 2000, nchains = 1, summary = TRUE, samplesAsCodaMCMC = TRUE)

cat("\n--- EXECUÇÃO CONCLUÍDA ---\n")
print("Sumário das amostras:")
print(mcmc.out_compiled$summary)

