source("_dataCaseStudy.R")
attach(data)
source("Covariaveis.R") 
source("ValoresIniciais.R")
source("dynamic_sampler_new.R")
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

n_clusters   <- 4
n_covariates <- as.integer(dim(x)[3])
n_regions    <- as.integer(dim(E)[1])
n_times      <- as.integer(dim(E)[2])

# A lista de constantes está correta. 
# a0 e b0 são usados pelo dynamic_sampler, não diretamente pelo nimbleCode.
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
  Y = Y,   # Certifique-se de que y_ini seja uma matriz de inteiros
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
samples <- runMCMC(Cmcmc, niter = 5000, nchains = nchains, inits = inits_list, samplesAsCodaMCMC = TRUE,summary = TRUE)

samples$summary

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
