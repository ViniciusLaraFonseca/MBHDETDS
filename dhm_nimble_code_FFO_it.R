# -------------------------------------------------------------------
# SCRIPT PRINCIPAL — TESTE (lambda_it, T=100) com FFO INTERNO
# (Filtro dentro do nimbleCode)
# -------------------------------------------------------------------

# --- PASSO 0: CARREGAR PACOTES E SCRIPTS ---
if (!require(nimble)) install.packages("nimble")
if (!require(coda)) install.packages("coda")
if (!require(ggplot2)) install.packages("ggplot2")
library(nimble)
library(coda)
library(ggplot2)

cat("--- Carregando dados e iniciais ---\n")
source("_dataCaseStudy.r") 
# source("ValoresIniciais_Lambda_it_T100.R") # Script 1

# --- PASSO 1: DEFINIÇÃO DO MODELO NIMBLE ---
cat("\n--- PASSO 1: Construindo modelo NIMBLE (FFO, lambda_it) ---\n")

code_ffo_it <- nimbleCode({
  
  # 1. Priors (Beta e Gamma)
  for (j in 1:p) {
    beta[j] ~ dnorm(mu_beta[j], sd = 5)
  }
  gamma[1] ~ dunif(a_unif[1], b_unif[1])
  for (j in 2:K) {
    gamma[j] ~ dunif(
      min = a_unif[j] * (1 - sum(gamma[1:(j - 1)])),
      max = b_unif[j] * (1 - sum(gamma[1:(j - 1)]))
    )
  }
  for (i in 1:n_regions) {
    epsilon[i] <- 1 - sum(h[i, 1:K] * gamma[1:K])
  }
  
  # 2. Componente Espacial g_it
  for (i in 1:n_regions) {
    for (t in 1:n_times) {
      log_g_it[i, t] <- inprod(beta[1:p], x[i, t, 1:p])
      g_it[i, t] <- E[i, t] * epsilon[i] * exp(log_g_it[i, t])
    }
  }
  
  # --- 3. FILTRO FORWARD (POR ÁREA) DENTRO DO MODELO ---
  for(i in 1:n_regions) {
    # Inicialização t=1
    at[i, 1] <- a0
    bt[i, 1] <- b0
    
    # Evolução t >= 2
    for (t in 2:(n_times + 1)) {
      # Preditiva
      att[i, t-1] <- w * at[i, t-1]
      btt[i, t-1] <- w * bt[i, t-1]
      
      # Atualização
      at[i, t] <- att[i, t-1] + Y[i, t-1]
      bt[i, t] <- btt[i, t-1] + g_it[i, t-1] 
    }
  }
  
  # --- 4. PRIORS PREDITIVAS + VEROSSIMILHANÇA ---
  for(i in 1:n_regions) {
    for (t in 1:n_times) {
      # Prior preditiva (forward-only)
      lambda[i, t] ~ dgamma(att[i, t], rate = btt[i, t])
      
      # Verossimilhança
      mu[i, t] <- lambda[i, t] * g_it[i, t]
      Y[i, t] ~ dpois(mu[i, t])
    }
  }
})

# --- PASSO 2: CONFIGURAÇÃO DO MCMC ---
cat("\n--- PASSO 2: Configurando MCMC (FFO, lambda_it) ---\n")

# Definir dimensões (precisa ser feito no R, não no nimbleCode)
dims_list = list(
  g_it      = c(constants_nimble$n_regions, constants_nimble$n_times),
  log_g_it  = c(constants_nimble$n_regions, constants_nimble$n_times),
  mu        = c(constants_nimble$n_regions, constants_nimble$n_times),
  lambda    = c(constants_nimble$n_regions, constants_nimble$n_times),
  at        = c(constants_nimble$n_regions, constants_nimble$n_times + 1),
  bt        = c(constants_nimble$n_regions, constants_nimble$n_times + 1),
  att       = c(constants_nimble$n_regions, constants_nimble$n_times),
  btt       = c(constants_nimble$n_regions, constants_nimble$n_times)
)

model_ffo_it <- nimbleModel(code_ffo_it,
                            constants = constants_nimble,
                            data      = data_nimble,
                            inits     = inits_list_nimble[[1]],
                            dimensions = dims_list, # Informa ao NIMBLE as dimensões dos buffers
                            check     = FALSE)

Cmodel_ffo_it <- compileNimble(model_ffo_it)

# Configurar MCMC (Sem sampler customizado)
conf_ffo_it <- configureMCMC(Cmodel_ffo_it, monitors = c("lambda", "beta", "gamma"))
conf_ffo_it$printSamplers()

Rmcmc_ffo_it <- buildMCMC(conf_ffo_it)
Cmcmc_ffo_it <- compileNimble(Rmcmc_ffo_it, project = Cmodel_ffo_it)

# --- PASSO 3: EXECUTAR MCMC ---
cat("\n--- Executando MCMC (lambda_it, FFO Interno, T=100) com 2 cadeias ---\n")
samples_ffo_it <- runMCMC(Cmcmc_ffo_it,
                          niter           = 10000,
                          nburnin         = 5000,
                          nchains         = 2,
                          inits           = inits_list_nimble,
                          samplesAsCodaMCMC = TRUE)

# --- PASSO 4: GERAR TRACEPLOTS ---
cat("\n--- Gerando traceplots (lambda_it, FFO Interno) ---\n")
# ... (use o mesmo código de plotagem) ...