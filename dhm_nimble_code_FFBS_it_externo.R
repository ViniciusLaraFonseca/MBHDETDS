# -------------------------------------------------------------------
# SCRIPT PRINCIPAL — TESTE (lambda_it, T=100) com FFBS EXTERNO
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
source("ffbs_sampler_Lambda_it.R")         # Script 2

# --- PASSO 1: DEFINIÇÃO DO MODELO NIMBLE ---
cat("\n--- PASSO 1: Construindo modelo NIMBLE ---\n")

code_ffbs_externo_it <- nimbleCode({
  
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
  
  # 3. Prior Fictício para lambda_it (MATRIZ 2D)
  for (i in 1:n_regions) {
    for (t in 1:n_times) {
      lambda[i, t] ~ dgamma(1, 1) # Será substituído pelo FFBS
    }
  }
  
  # 4. Verossimilhança
  for (i in 1:n_regions) {
    for (t in 1:n_times) {
      mu[i, t] <- lambda[i, t] * g_it[i, t] # Usa lambda[i,t]
      Y[i, t] ~ dpois(mu[i, t])
    }
  }
})

# --- PASSO 2: CONFIGURAÇÃO DO MCMC ---
cat("\n--- PASSO 2: Configurando MCMC ---\n")

model_ext_it <- nimbleModel(code_ffbs_externo_it,
                            constants = constants_nimble,
                            data      = data_nimble,
                            inits     = inits_list_nimble[[1]],
                            check     = FALSE)

Cmodel_ext_it <- compileNimble(model_ext_it)

conf_ext_it <- configureMCMC(model = Cmodel_ext_it,
                             monitors = c("lambda", "beta", "gamma"))

cat("Substituindo amostrador de lambda pelo FFBS (por área)...\n")
conf_ext_it$removeSamplers("lambda")
conf_ext_it$addSampler(target = "lambda",
                       type = "ffbs_sampler_Lambda_it", # O novo sampler 2D
                       control = list(
                         n_regions = constants_nimble$n_regions,
                         n_times   = constants_nimble$n_times,
                         p         = constants_nimble$p,
                         w         = constants_nimble$w,
                         a0        = constants_nimble$a0,
                         b0        = constants_nimble$b0
                       ))
conf_ext_it$printSamplers()

Rmcmc_ext_it <- buildMCMC(conf_ext_it)
Cmcmc_ext_it <- compileNimble(Rmcmc_ext_it, project = Cmodel_ext_it)

# --- PASSO 3: EXECUTAR MCMC ---
cat("\n--- Executando MCMC (lambda_it, FFBS Externo, T=100) com 2 cadeias ---\n")
samples_ext_it <- runMCMC(Cmcmc_ext_it,
                          niter           = 10000, # Reduzido para teste (T=100 é longo)
                          nburnin         = 5000,
                          nchains         = 2,
                          inits           = inits_list_nimble,
                          samplesAsCodaMCMC = TRUE)

# --- PASSO 4: GERAR TRACEPLOTS ---
cat("\n--- Gerando traceplots (lambda_it, FFBS Externo) ---\n")
# ... (use o mesmo código de plotagem do script anterior) ...