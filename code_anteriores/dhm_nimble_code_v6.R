code <- nimbleCode({
  # Priors para coeficientes de regressão
  for (j in 1:p) {
    beta[j] ~ dnorm(mu_beta[j], sd = 5)
  }
  
  # Priors para parâmetros de subnotificação
  gamma[1] ~ dunif(a_unif[1], b_unif[1])
  for (j in 2:K) {
    gamma[j] ~ dunif(
      min = a_unif[j] * (1 - sum(gamma[1:(j - 1)])),
      max = b_unif[j] * (1 - sum(gamma[1:(j - 1)]))
    )
  }
  
  # Estrutura para probabilidade de notificação
  for (i in 1:A) {
    epsilon[i] <- 1 - sum(h[i, 1:K] * gamma[1:K])
  }
  
  # Modelo para estados latentes e observações
  for (i in 1:A) {
    # Inicialização do processo dinâmico
    lambda[i, 1] ~ dgamma(a0, b0)
    log_theta[i, 1] <- log(lambda[i, 1]) + inprod(beta[1:p], x[i, 1, 1:p])
    theta[i, 1] <- exp(log_theta[i, 1])
    
    # Equação de observação
    mu[i,1] <- E[i, 1] * epsilon[i] * theta[i, 1]
    Y[i, 1] ~ dpois(mu[i,1])
    
    # Evolução temporal
    for (t in 2:T) {
      # Equação de estado para lambda
      lambda[i, t] ~ dgamma(w * a_prev[i, t-1], w * b_prev[i, t-1])
      
      # Risco relativo
      log_theta[i, t] <- log(lambda[i, t]) + inprod(beta[1:p], x[i, t, 1:p])
      theta[i, t] <- exp(log_theta[i, t])
      
      # Equação de observação
      mu[i,t] <- E[i, t] * epsilon[i] * theta[i, t]
      Y[i, t] ~ dpois(mu[i,t])
      
      # Atualização dos parâmetros
      a_prev[i, t] <- w*a_prev[i, t-1] + Y[i, t]
      b_prev[i, t] <- w*b_prev[i, t-1] + E[i, t] * epsilon[i] * exp(inprod(beta[1:p], x[i,t,1:p]))
    }
    
    # Inicialização dos parâmetros de filtragem
    a_prev[i, 1] <- a0 + Y[i, 1]
    b_prev[i, 1] <- b0 + E[i, 1] * epsilon[i] * exp(inprod(beta[1:p], x[i,1,1:p]))
  }
})

n_clusters = 4
n_covariates = dim(x)[3]
n_regions = dim(E)[1]  # Número de regiões
n_times = dim(E)[2] 
h_matrix = hAI
Y_matrix = y_ini
E_matrix = E
x_array = x
constants <- list(
  mu_beta = beta_ini,
  w=0.9,
  A = n_regions,
  T=n_times,
  a_unif = a_unif_ini,
  b_unif = b_unif_ini,
  n_regions = n_regions,
  n_times = n_times,
  p = n_covariates,
  K = n_clusters,
  a0 = 1,    # hiperparâmetros iniciais
  b0 = 1,
  h = h_matrix  # matriz de indicadores de cluster
)

data <- list(
  Y = round(Y_matrix),      # dados observados
  E = E_matrix,      # valores esperados
  x = x_array        # array de covariáveis [A, T, p]
)

# Inicialização
inits <- list(
  beta = rnorm(constants$p, 0, 1),
  gamma = gamma_ini,
  lambda = matrix(lambda0, constants$n_regions, constants$n_times)
)

# Construir modelo
model <- nimbleModel(code, constants = constants, data = data, inits = inits)

# Configurar MCMC
conf <- configureMCMC(model, monitors = c("beta", "gamma", "lambda", "theta"))

# Adicionar amostrador customizado para lambda
conf$removeSamplers("lambda","beta","gamma")
conf$addSampler(target = "beta", type = "RW_block")
conf$addSampler(target = "gamma", type = "RW_block")
conf$addSampler(target = "lambda", type = "dynamic_sampler", 
                control = list(n_regions = constants$n_regions,n_times = constants$n_times,w =constants$w))
conf$samplerConfs
# Compilar e executar
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(model)
Cmcmc <- compileNimble(Rmcmc, project = model)

nchains = 2
samples <- runMCMC(Cmcmc,nchains =nchains ,niter = 25000)

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

calcNodes <- model$getDependencies("lambda")    
nodes <- calcNodes
# pega log-prob por nó
lp_by_node <- sapply(nodes, function(n) model$getLogProb(n))
# ordena e mostra os 20 piores
ord <- order(lp_by_node)
data.frame(node = nodes[ord[1:20]], logProb = lp_by_node[ord[1:20]])

range(values(Cmodel, 'epsilon'), na.rm=TRUE)
range(values(Cmodel, 'theta'), na.rm=TRUE)
range(values(Cmodel, 'mu'), na.rm=TRUE)
range(values(Cmodel, 'lambda'), na.rm=TRUE)
range(values(Cmodel, 'a_prev'), na.rm=TRUE)
range(values(Cmodel, 'b_prev'), na.rm=TRUE)

i <- 74
T <- constants$T
A <- constants$A
idx <- ((i - 1) * T + 1):(i * T)
values(Cmodel, 'Y')[idx]
