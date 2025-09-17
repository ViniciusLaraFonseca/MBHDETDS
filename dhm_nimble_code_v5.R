K = 4

# Estimando os coeficientes iniciais usando um modelo GLM (Modelo Linear Generalizado)
# O modelo utiliza a variável resposta Y2009_2011 e preditores selecionados da matriz x.
n_regions = dim(E)[1]  # Número de regiões
n_times = dim(E)[2]     # Número de períodos de tempo


glm <- glm(Y2009_2011 ~ 1 + x[, 10, 1] , family = poisson(link = "log"))
summary(glm)
# Extraindo o coeficiente intercepto (lambda0) do modelo, para calcular a0 e b0
lambda0 <- glm$coefficients[1]

# Extraindo os coeficientes beta iniciais (para as variáveis preditoras)
beta_ini <- glm$coefficients[2]
# Definindo valores iniciais para os parâmetros a e b
a0 =18.7321  # Parâmetro a inicial
b0 = 4.32795  # Parâmetro b inicial

# Definindo valores iniciais para os parâmetros gamma
gamma_ini = c(0.05, 0.10, 0.10, 0.15)

# Inicializando os parâmetros a e b para uma distribuição uniforme
a_unif_ini = c(0.00, rep(0, K - 1))  # Parâmetros a uniformes
b_unif_ini = c(0.1, rep(1, K - 1))   # Parâmetros b uniformes

# Definindo valores iniciais para os parâmetros gamma
gamma_ini = c(0.05, 0.10, 0.10, 0.15)

# Inicializando os parâmetros a e b para uma distribuição uniforme
a_unif_ini = c(0.00, rep(0, K - 1))  # Parâmetros a uniformes
b_unif_ini = c(0.1, rep(1, K - 1))   # Parâmetros b uniformes

# Estimando epsilon com base em um fator hAI e os parâmetros gamma
epsilon_ini = 1 - hAI %*% gamma_ini
# Expandindo epsilon para uma matriz com 75 regiões e 23 períodos
epsilon_ini = matrix(rep(epsilon_ini, 23), nrow = 75, ncol = 23)

# Estimando lambda inicial utilizando a distribuição gamma
lambda_ini1 = matrix(rgamma(1725, shape = a0, rate = b0), ncol = 23, nrow = 75)
lambda_ini2 = matrix(rep(1,1725), ncol = 23, nrow = 75)
lambda_ini3 = matrix(rep(2,1725), ncol = 23, nrow = 75)
lambda_ini4 = matrix(rep(5,1725), ncol = 23, nrow = 75)
lambda_ini5 = matrix(rep(10,1725), ncol = 23, nrow = 75)
#estimando o theta com os betas iniciais e o x
theta_ini = lambda_ini1 * exp(beta_ini[1] * x[,, 1] )

# Estimando Y com base nos parâmetros lambda, E e epsilon, aplicando a função exponencial nos coeficientes beta
y_ini =  round(E * epsilon_ini * theta_ini)

#estimando o theta com os betas iniciais e o x


# Definindo um parâmetro w para ser utilizado em cálculos posteriores
w=0.7
code <- nimbleCode({
  # Priors for beta coefficients (covariates in theta)
  
  beta ~ dnorm(mu_beta, sd = sd_beta)
  
  
  # Priors for gamma (used in epsilon computation)
  gamma[1] ~ dunif(a_unif[1], b_unif[1])
  for (j in 2:K) {
    gamma[j] ~ dunif(
      min = a_unif[j] * (1 - sum(gamma[1:(j - 1)])),
      max = b_unif[j] * (1 - sum(gamma[1:(j - 1)]))
    )
  }
  for (t in 1:n_times){
    for (i in 1:n_regions) {
      
      # Efeito das covariáveis qualitativas em epsilon
      epsilon[i, t] <- 1 - inprod(hAI[i, 1:K], gamma[1:K])
      
      # Covariado contínuo em theta
      theta[i, t] <- lambda[i, t]*exp((beta* x[i, t]))
      
      # Esperança da Poisson
      mu[i, t] <-  E[i, t] * epsilon[i, t] * theta[i, t]
      
      # Likelihood
      Y[i, t] ~ dpois(mu[i, t])
    }
  }
  for (i in 1:n_regions) {
    for (t in 1:n_times) {
      lambda[i, t] ~ dgamma(1, 1)  # dummy prior
    }
  }
}
)

modelConstants = list(K= K,
                      n_regions = dim(x)[1],
                      n_times = dim(x)[2],
                      n_covariates_theta = 1, 
                      sd_beta = 5,
                      a_unif = a_unif_ini,
                      b_unif = b_unif_ini,
                      a0 = a0,
                      b0 = b0,
                      mu_beta = beta_ini,
                      count =y_ini,
                      w=w
)
modelData = list(Y = y_ini,
                 x = x[,,1],
                 hAI = hAI,
                 E = E)
modelInits1 = list(lambda = lambda_ini1,
                   gamma = gamma_ini,
                   beta = beta_ini)
modelInits2 = list(lambda = lambda_ini2,
                   gamma = c(0.10,0.15,0.4,0.6),
                   beta = -0.5)
modelInits3 = list(lambda = lambda_ini3,
                   gamma = c(0.14,0.10,0.15,0.5),
                   beta = 1)
modelInits4 = list(lambda = lambda_ini4,
                   gamma = c(0.10,0.10,0.10,0.10),
                   beta = -1)

modelInits <- c(modelInits1,modelInits2,modelInits3,modelInits4)
#Create the model
#modelInits = c(modelInits1,modelInits2)
model = nimbleModel(code,constants = modelConstants,data = modelData,inits = modelInits)
#compilando o modelo
cModel = compileNimble(model)

#configurando os targets para o modelo 
modelConf = configureMCMC(model,monitors=c('beta','gamma',
                                           'mu',"lambda",'epsilon',
                                           'theta'), useConjugacy = TRUE,enableWAIC = F)
#removendo e adicionando o sampler
modelConf$removeSamplers("lambda","gamma")
modelConf$addSampler(target = "lambda",type = "dynamic_sampler",control = list(
  E = modelData$E,
  x = modelData$x,
  w = modelConstants$w,
  a0 = modelConstants$a0,
  b0 = modelConstants$b0,
  n_regions = modelConstants$n_regions,
  n_times = modelConstants$n_times,
  count = modelConstants$count
))
modelConf$addSampler(
  target = paste0("gamma[", 1:K, "]"),
  type   = "RW_block",
  control = list(
    adaptScale = TRUE,   # adapta a escala automaticamente
    adaptInterval = 200, # intervalo de adaptação
    scale = 0.1  # tamanho inicial dos passos
  )
)
modelConf$samplerConfs
modelMCMC = buildMCMC(modelConf)

cMCMC <- compileNimble(modelMCMC, project = model, resetFunctions = TRUE)
#definindo o número de cadeias
nchains = 1
#rodando as cadeias
mcmc.out = runMCMC(cMCMC,
                   nchains =nchains, niter = 50000,
                   summary = TRUE, WAIC = FALSE)
#alocando as amostras
samples <- mcmc.out$samples
mcmc_list <- mcmc.list(lapply(1:nchains,function(chain){
  mcmc(samples[[chain]])
}))


beta_samples<- lapply(mcmc_list,function(chain){
  as.mcmc(chain)[,"beta"]
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
lambda1_10_samples <- lapply(mcmc_list, function(chain) {
  as.mcmc(chain)[, "lambda[1, 10]"]
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

par(mfrow = c(2,2))
for(i in 1:4){
  traceplot(lambda1_1_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
  abline(h = lambda_ini1[1,1],col = "red")
}
par(mfrow = c(2,2))
for(i in 1:4){
  traceplot(lambda1_2_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
  abline(h = lambda_ini1[1,2],col = "red")
}
par(mfrow = c(2,2))
for(i in 1:4){
  traceplot(lambda1_3_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
  abline(h = lambda_ini1[1,3],col = "red")
}
par(mfrow = c(2,2))
for(i in 1:4){
  traceplot(lambda1_10_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
  abline(h = lambda_ini1[1,10],col = "red")
}
par(mfrow = c(2,2))
for(i in 1:4){
  traceplot(lambda1_23_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
  abline(h = lambda_ini1[1,23],col = "red")
}
par(mfrow = c(2,2))
for(i in 1:4){
  traceplot(gamma1_samples[i],main=paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
  abline(h = gamma_ini[1],col = "red")
}
par(mfrow = c(2,2))
for(i in 1:4){
  traceplot(gamma2_samples[i],main=paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
  abline(h = gamma_ini[2],col = "red")
}
par(mfrow = c(2,2))
for(i in 1:4){
  traceplot(gamma3_samples[i],main=paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
  abline(h = gamma_ini[3],col = "red")
}
par(mfrow = c(2,2))
for(i in 1:4){
  traceplot(gamma4_samples[i],main=paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
  abline(h = gamma_ini[4],col = "red")
}
par(mfrow = c(2,2))
for(i in 1:nchains){
  traceplot(beta_samples[i],main = paste("Cadeia",i),xlab = "Iterações",ylab = "Valores")
  abline(h = beta_ini[1],col = "red")
}
