# -------------------------------
# 1️⃣ Código do modelo
# -------------------------------
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
      log_theta[i, t] <- log(lambda[i, t]) + inprod(beta[1:p], x[i, t, 1:p])
      theta[i, t] <- exp(log_theta[i, t])
      mu[i,t] <- E[i, t] * epsilon[i] * theta[i, t]
      Y[i, t] ~ dpois(mu[i,t])
    }
  }
})

# -------------------------------
# 2️⃣ Constantes e dados
# -------------------------------
# Constantes e Dados
# -------------------------------
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

# -------------------------------
# 5️⃣ Configurar MCMC
# -------------------------------

model <- nimbleModel(code, constants = constants, data = data, inits = inits1)
conf <- configureMCMC(model, monitors = c("beta", "gamma", "lambda", "theta"))
conf$removeSamplers("lambda")
conf$addSampler(target = "lambda", type = dynamic_sampler,
                control = list(w = 0.9,
                               a0 = constants$a0,
                               b0 = constants$b0))
print(conf$samplerConfs)

# -------------------------------
# 6️⃣ Compilar modelo e MCMC
# -------------------------------

Cmodel <- compileNimble(model,resetFunctions = TRUE)
Rmcmc  <- buildMCMC(conf)
Cmcmc  <- compileNimble(Rmcmc, project = Cmodel, showCompilerOutput = TRUE)

# -------------------------------
#compilar modelo e MCMC sem o dynamic_sampler
# -------------------------------

conf_test <- configureMCMC(model, monitors = c("beta", "gamma", "lambda", "theta"))
Rmcmc_test <- buildMCMC(conf_test)
Cmcmc  <- compileNimble(Rmcmc_test, project = Cmodel, showCompilerOutput = TRUE,resetFunctions=TRUE)

# -------------------------------
# 7️⃣ Rodar MCMC com múltiplas cadeias
# -------------------------------
nchains = length(inits_list)
samples <- runMCMC(Cmcmc, niter = 5000, nchains = nchains, inits = inits_list, samplesAsCodaMCMC = TRUE)


mcmc_list <- mcmc.list(lapply(1:nchains,function(chain){
  mcmc(samples[[chain]])
}))


beta_1_samples<- lapply(mcmc_list,function(chain){
  as.mcmc(chain)[,"beta[1]"]
})
beta_2_samples<- lapply(mcmc_list,function(chain){
  as.mcmc(chain)[,"beta[2]"]
})
beta_3_samples<- lapply(mcmc_list,function(chain){
  as.mcmc(chain)[,"beta[3]"]
})
lambda1_1_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 1]"]
})
lambda1_2_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 2]"]
})
lambda1_3_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 3]"]
})
lambda1_4_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 4]"]
})
lambda1_5_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 5]"]
})
lambda1_10_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 10]"]
})
lambda1_11_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 11]"]
})
lambda1_12_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 12]"]
})
lambda1_13_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 13]"]
})
lambda1_14_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 14]"]
})
lambda1_15_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 15]"]
})
lambda1_16_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 16]"]
})
lambda1_17_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 17]"]
})
lambda1_18_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 18]"]
})
lambda1_19_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 19]"]
})
lambda1_20_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 20]"]
})
lambda1_21_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 21]"]
})
lambda1_22_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 22]"]
})
lambda1_23_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 23]"]
})

lambda2_23_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[2, 23]"]
})

lambda2_1_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[2, 1]"]
})
gamma1_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "gamma[1]"]})
gamma2_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "gamma[2]"]})
gamma3_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "gamma[3]"]})
gamma4_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "gamma[4]"]})

#traceplot lambda
par(mfrow = c(2,1))
for(i in 1:nchains){
  traceplot(lambda1_1_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
  
}

for(i in 1:nchains){
  traceplot(lambda1_2_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
  
}

for(i in 1:nchains){
  traceplot(lambda1_3_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
}

for(i in 1:nchains){
  traceplot(lambda1_10_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
}
for(i in 1:nchains){
  traceplot(lambda1_11_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
}
for(i in 1:nchains){
  traceplot(lambda1_12_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
}
for(i in 1:nchains){
  traceplot(lambda1_13_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
}
for(i in 1:nchains){
  traceplot(lambda1_14_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
}
for(i in 1:nchains){
  traceplot(lambda1_15_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
}
for(i in 1:nchains){
  traceplot(lambda1_16_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
}
for(i in 1:nchains){
  traceplot(lambda1_17_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
}
for(i in 1:nchains){
  traceplot(lambda1_18_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
}
for(i in 1:nchains){
  traceplot(lambda1_19_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
}
for(i in 1:nchains){
  traceplot(lambda1_20_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
}
for(i in 1:nchains){
  traceplot(lambda1_21_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
}
for(i in 1:nchains){
  traceplot(lambda1_22_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
}
for(i in 1:nchains){
  traceplot(lambda1_23_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
}
for(i in 1:nchains){
  traceplot(lambda2_23_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
}
for(i in 1:nchains){
  traceplot(lambda2_1_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
}
#traceplot betas
for(i in 1:nchains){
  traceplot(beta_1_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")}

for(i in 1:nchains){
  traceplot(beta_2_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")}

for(i in 1:nchains){
  traceplot(beta_3_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")}

#traceplot gammas

for(i in 1:nchains){
  traceplot(gamma1_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")}
for(i in 1:nchains){
  traceplot(gamma2_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")}
for(i in 1:nchains){
  traceplot(gamma3_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")}
for(i in 1:nchains){
  traceplot(gamma4_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")}

