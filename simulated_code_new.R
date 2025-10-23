source("_dataCaseStudy.R")
code <- nimbleCode({
  # Priors para os coeficientes de regressão
  for (j in 1:p) {
    beta[j] ~ dnorm(mu_beta[j], sd = 5)
  }
  
  # Priors para os parâmetros de subnotificação
  gamma[1] ~ dunif(a_unif[1], b_unif[1])
  for (j in 2:K) {
    gamma[j] ~ dunif(
      min = a_unif[j] * (1 - sum(gamma[1:(j - 1)])),
      max = b_unif[j] * (1 - sum(gamma[1:(j - 1)]))
    )
  }
  
  # Probabilidade de notificação
  for (i in 1:n_regions) {
    epsilon[i] <- 1 - sum(h[i, 1:K] * gamma[1:K])
  }
  
  # Verossimilhança para as observações
  for (i in 1:n_regions) {
    for(t in 1:n_times){  
      # Prior Fictício: Apenas para que NIMBLE reconheça lambda como estocástico
      lambda[i, t] ~ dgamma(1, 1)
      log_theta[i, t] <- log(lambda_suv[i, t]) + inprod(beta[1:p], x[i, t, 1:p])
      theta[i, t] <- exp(log_theta[i, t])
      mu[i,t] <- E[i, t] * epsilon[i] * theta[i, t]
      Y[i, t] ~ dpois(mu[i,t])
    }
    lambda_suv[i,n_times] <- lambda[i,n_times]
    for(t in 1:(n_times-1)) {
      lambda_suv[i, n_times-t] <- lambda[i,n_times-t] + w * lambda[i,n_times-t+1]
    }
  }
})
set.seed(123)
#apesar de ser uma simulação vamos manter os dados originais
n_clusters   <- 4
n_covariates <- 3
n_regions    <- 75
n_times      <- 23

#assumir a forma mais simples de lambda0
a0=1
b0=1

E_at <- matrix(runif(n_regions * (n_times+1), 100, 200), nrow = n_regions)
E <- matrix(E_at[,2:(n_times+1)],nrow=n_regions)
x_at <- array(rnorm(n_regions * (n_times+1) * p), dim = c(n_regions, (n_times+1), p))
x <- array(x_at[,2:(n_times+1),],dim = c(n_regions, (n_times+1), p))
hAI <- data$hAI
epsilon <- 1 - as.numeric(hAI %*% gamma)

# Salvar valores reais de lambda para comparação

lambda0_real <- rgamma(n_times, shape = a0, rate = b0)

#EFEITO DINAMICO 
lambda0 <- lambda0_real
mu_true0<- E[1,1] * epsilon * exp(c(x[1,1,] %*% beta_initial_values)) * lambda0
y_ini0 <- matrix(rpois(1, lambda = mu_true0), nrow = n_regions)
lambda_at_fwd <- matrix(NA, ncol = n_times, nrow = n_regions)
y_at <- matrix(NA, nrow = n_regions, ncol = n_times + 1)
y_at[, 1] <- y_ini0

att_ini <- matrix(NA, nrow = n_regions, ncol = n_times)
btt_ini <- matrix(NA, nrow = n_regions, ncol = n_times)
at_ini  <- matrix(NA, nrow = n_regions, ncol = n_times + 1)
bt_ini  <- matrix(NA, nrow = n_regions, ncol = n_times + 1)

# Simulação Forward
for(i in 1:n_regions) {
  at_ini[i, 1] <- a0
  bt_ini[i, 1] <- b0
  for(t in 2:(n_times + 1)) {
    att_ini[i, t-1] <- w * at_ini[i, t-1]
    btt_ini[i, t-1] <- w * bt_ini[i, t-1]
    at_ini[i, t]    <- att_ini[i, t-1] + y_at[i, t-1]
    
    prod_val <- sum(x_at[i, t, ] * beta_initial_values)
    bt_ini[i, t]    <- btt_ini[i, t-1] + E_at[i, t] * epsilon[i] * exp(prod_val)
    
    lambda_at_fwd[i,t-1] <- rgamma(1, shape = at_ini[i, t], rate = bt_ini[i, t])
    
    mu_at <- lambda_at_fwd[i, t-1] * exp(prod_val) * epsilon[i] * E_at[i, t]
    y_at[i, t] <- rpois(1, mu_at)
  }
}
y_ini <- y_at[,2:21, drop = FALSE]

dynamic_sampler <- nimbleFunction(
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    # Extrair dimensões e parâmetros
    dims_Y <- dim(model$Y)
    n_regions <- dims_Y[1]
    n_times   <- dims_Y[2]
    p         <- dim(model$x)[3]
    w         <- control$w
    a0        <- control$a0
    b0        <- control$b0
    
    # Buffers para armazenar os parâmetros da filtragem
    at_buf  <- nimMatrix(nrow = n_regions, ncol = n_times + 1, init = 0, type = 'double')
    bt_buf  <- nimMatrix(nrow = n_regions, ncol = n_times + 1, init = 0, type = 'double')
    
    # Nós a serem calculados após a amostragem
    calcNodes   <- model$getDependencies(target, self = FALSE)
    targetNodes <- model$expandNodeNames(target)
  },
  
  run = function() {
    # --- Declaração explícita de tipos ---
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
    
    # --- Forward Filtering ---
    for(i in 1:n_regions) {
      at_buf[i, 1] <<- a0
      bt_buf[i, 1] <<- b0
      
      for(t in 1:n_times) {
        att_t <- w * at_buf[i, t]
        btt_t <- w * bt_buf[i, t]
        
        at_buf[i, t+1] <<- att_t + model$Y[i, t]
        
        prod_val <- 0
        for(k in 1:p) {
          prod_val <- prod_val + model$x[i, t, k] * model$beta[k]
        }
        
        bt_buf[i, t+1] <<- btt_t + model$E[i, t] * model$epsilon[i] * exp(prod_val)
      }
    }
    
    # --- Backward Sampling ---
    for(i in 1:n_regions) {
      shape_tmp <- at_buf[i, n_times + 1]
      rate_tmp  <- bt_buf[i, n_times + 1]
      model$lambda[i, n_times] <<- rgamma(1, shape = shape_tmp, rate = rate_tmp)
      
      for(t in 1:(n_times-1)) {
        lambda_futuro <- model$lambda[i, (n_times-t+1)]
        
        shape_tmp <- (1 - w) * at_buf[i, (n_times-t+1)] 
        rate_tmp  <- bt_buf[i, t]
        
        model$lambda[i,(n_times-t)] <<- rgamma(1, shape = shape_tmp, rate = rate_tmp)
        
        model$lambda_suv[i, (n_times-t)] <<- model$lambda[i,(n_times-t)] + w * lambda_futuro
      }
    }
    
    # Atualizar os nós dependentes e a log-probabilidade
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = TRUE)
  },
  
  methods = list(
    reset = function() {}
  )
)

constants <- list(
  mu_beta   = beta_ini,
  w         = 0.9,
  a_unif    = a_unif_ini,
  b_unif    = b_unif_ini,
  n_regions = n_regions,
  n_times   = n_times,
  p         = n_covariates,
  K         = n_clusters,
  a0        = a0,
  b0        = b0,
  h         = hAI
)
print(constants)

# A lista de dados foi simplificada.
# A variável 'count' é redundante, pois o modelo e o amostrador usam 'Y'.
data <- list(
  Y = y_ini,   # Certifique-se de que y_ini seja uma matriz de inteiros
  E = E,
  x = x
)
print(data)

# -------------------------------
# 3️⃣ Inicializações para 2 cadeias
# -------------------------------
# A estrutura das inicializações está correta.
# É importante fornecer valores iniciais para o estado latente 'lambda'.
inits1 <- list(
  beta   = rnorm(constants$p, beta_ini, 1),
  gamma  = gamma_ini,
  lambda = lambda_ini   # Assumindo que 'lambda_ini' já é uma matriz
)
print(inits1)

inits2 <- list(
  beta   = rnorm(constants$p, 0, 0.5),
  gamma  = gamma_ini / 2,
  lambda = matrix(1, n_regions, n_times)
)

inits_list <- list(inits1, inits2)

model <- nimbleModel(code, constants = constants, data = data, inits = inits1)
conf <- configureMCMC(model, monitors = c("beta", "gamma", "lambda","lambda_suv", "theta"))
conf$removeSamplers("lambda")
conf$addSampler(target = "lambda", type = dynamic_sampler,
                control = list(w = 0.9,
                               a0 = constants$a0,
                               b0 = constants$b0))
print(conf$samplerConfs)

Cmodel <- compileNimble(model,resetFunctions = TRUE)
Rmcmc  <- buildMCMC(conf)
Cmcmc  <- compileNimble(Rmcmc, project = Cmodel, showCompilerOutput = FALSE)

nchains = length(inits_list)
samples <- runMCMC(Cmcmc, niter = 5000, nchains = nchains, inits = inits_list, samplesAsCodaMCMC = TRUE)


mcmc_list <- mcmc.list(lapply(1:nchains,function(chain){
  mcmc(samples[[chain]])
}))

for (i in 1:constants$p) {
  png(paste0("traceplot_beta_", i, ".png"))
  coda::traceplot(mcmc_list[, paste0("beta[", i, "]")], main = paste("Traceplot for beta[", i, "]"))
  dev.off()
}

# Parâmetros gamma
for (i in 1:constants$K) {
  png(paste0("traceplot_gamma_", i, ".png"))
  coda::traceplot(mcmc_list[, paste0("gamma[", i, "]")], main = paste("Traceplot for gamma[", i, "]"))
  dev.off()
}

# Parâmetros lambda selecionados (baseado no seu relatório de diagnóstico)
selected_lambdas <- c("lambda[1, 1]", "lambda[1, 19]","lambda[1, 23]")
for (lam in selected_lambdas) {
  filename <- paste0("traceplot_", gsub("\\[|, |\\]", "_", lam), ".png")
  png(filename)
  coda::traceplot(mcmc_list[, lam], main = paste("Traceplot for", lam))
  dev.off()
}

selected_lambdas_suv <- c("lambda_suv[1, 1]", "lambda_suv[1, 19]","lambda_suv[1, 23]")
for (lam in selected_lambdas_suv) {
  filename <- paste0("traceplot_", gsub("\\[|, |\\]", "_", lam), ".png")
  png(filename)
  coda::traceplot(mcmc_list[, lam], main = paste("Traceplot for", lam))
  dev.off()
}