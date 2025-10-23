if (!require(nimble)) {
  install.packages("nimble")
}
if (!require(coda)) {
  install.packages("coda")
}
library(nimble)
library(coda)
# -------------------------------
# 1️⃣ Código do modelo com lambda transformado
# -------------------------------
code_transformed <- nimbleCode({
  # Priors para os coeficientes de regressão (beta) - Sem alterações
  for (j in 1:p) {
    beta[j] ~ dnorm(mu_beta[j], sd = 5)
  }
  
  # Priors para os parâmetros de subnotificação (gamma) - Sem alterações
  gamma[1] ~ dunif(a_unif[1], b_unif[1])
  for (j in 2:K) {
    gamma[j] ~ dunif(
      min = a_unif[j] * (1 - sum(gamma[1:(j - 1)])),
      max = b_unif[j] * (1 - sum(gamma[1:(j - 1)]))
    )
  }
  
  # Probabilidade de notificação (epsilon) - Sem alterações
  for (i in 1:n_regions) {
    epsilon[i] <- 1 - sum(h[i, 1:K] * gamma[1:K])
  }
  
  # Verossimilhança com o nó lambda_tilde explícito
  for (i in 1:n_regions) {
    for(t in 1:n_times){  
      # Prior fictício para lambda (o FFBS fará a amostragem real)
      lambda[i, t] ~ dgamma(1, 1)
      
      # --- TRANSFORMAÇÃO ---
      # 1. Define o componente g_t da transformação na escala log
      log_g[i,t] <- inprod(beta[1:p], x[i, t, 1:p])
      # 2. Cria o nó determinístico lambda_tilde
      lambda_tilde[i,t] <- lambda[i,t] * exp(log_g[i,t])
      
      # 3. A média da Poisson agora depende de lambda_tilde
      mu[i,t] <- E[i, t] * epsilon[i] * lambda_tilde[i, t]
      
      # 4. A verossimilhança
      Y[i, t] ~ dpois(mu[i,t])
    }
  }
})

# Carregar scripts necessários
source("_dataCaseStudy.R")
attach(data)
source("Covariaveis.R") 
source("ValoresIniciais.R")
source("dynamic_sampler.R") # O mesmo amostrador é utilizado

# -------------------------------
# 2️⃣ Constantes e Dados (sem alterações)
# -------------------------------
n_clusters   <- 4
n_covariates <- as.integer(dim(x)[3])
n_regions    <- as.integer(dim(E)[1])
n_times      <- as.integer(dim(E)[2])
nchains <- 2
constants <- list(
  mu_beta   = beta_ini,
  w         = 0.9,
  a_unif    = a_unif_ini,
  b_unif    = b_unif_ini,
  n_regions = n_regions,
  n_times   = n_times,
  p         = n_covariates,
  K         = n_clusters,
  a0        = 1,
  b0        = 1,
  h         = hAI
)

data <- list(Y = y_ini, E = E, x = x)

# -------------------------------
# 3️⃣ Inicializações (sem alterações)
# -------------------------------
inits1 <- list(
  beta   = rnorm(constants$p, beta_ini, 1),
  gamma  = gamma_ini,
  lambda = lambda_ini
)
inits2 <- list(
  beta   = rnorm(constants$p, 0, 0.5),
  gamma  = gamma_ini / 2,
  lambda = matrix(1, n_regions, n_times)
)
inits_list <- list(inits1, inits2)

# -------------------------------
# 4️⃣ Construção e Configuração do MCMC
# -------------------------------
cat("Construindo o modelo com lambda transformado...\n")
model_transformed <- nimbleModel(code_transformed, constants = constants, data = data, inits = inits1)

cat("Configurando o MCMC...\n")
conf_transformed <- configureMCMC(model_transformed, monitors = c("beta", "gamma", "lambda", "lambda_tilde"))
conf_transformed$removeSamplers("lambda")

# Adiciona o mesmo amostrador personalizado, pois a matemática para lambda não muda
conf_transformed$addSampler(target = "lambda", type = dynamic_sampler,
                            control = list(w = constants$w,
                                           a0 = constants$a0,
                                           b0 = constants$b0))

print(conf_transformed$samplerConfs)

# -------------------------------
# 5️⃣ Compilação e Execução
# -------------------------------
cat("Compilando o modelo e o MCMC...\n")
Cmodel_transformed <- compileNimble(model_transformed, resetFunctions = TRUE)
Rmcmc_transformed  <- buildMCMC(conf_transformed)
Cmcmc_transformed  <- compileNimble(Rmcmc_transformed, project = Cmodel_transformed, showCompilerOutput = TRUE)

cat("Executando o MCMC...\n")
samples_transformed <- runMCMC(Cmcmc_transformed, niter = 5000, nchains = 2, inits = inits_list, samplesAsCodaMCMC = TRUE)
mcmc_list <- coda::mcmc.list(lapply(1:nchains, function(chain) {
  coda::mcmc(samples_transformed[[chain]])
}))
cat("Gerando e salvando os traceplots...\n")

# Parâmetros beta
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
