code <- nimbleCode({
  
  # Priors for beta coefficients (covariates in theta)
  for (l in 1:n_covariates_theta) {
    beta[l] ~ dnorm(mu_beta[l], sd = sd_beta)
  }
  
  # Loop over time
  for (t in 1:n_times) {
    
    # Priors for gamma (used in epsilon computation)
    gamma[1, t] ~ dunif(a_unif[1], b_unif[1])
    for (j in 2:K) {
      gamma[j, t] ~ dunif(
        min = a_unif[j] * (1 - sum(gamma[1:(j - 1), t])),
        max = b_unif[j] * (1 - sum(gamma[1:(j - 1), t]))
      )
    }
    
    # Loop over regions
    for (i in 1:n_regions) {
      
      # Efeito das covariáveis qualitativas em epsilon
      epsilon[i, t] <- 1 - inprod(hAI[i, 1:K], gamma[1:K, t])
      
      # Covariado contínuo em theta
      theta[i, t] <- lambda[i, t]*exp(inprod(beta[1:n_covariates_theta], x[i, t, 1:n_covariates_theta]))
      
      # Esperança da Poisson
      mu[i, t] <-  E[i, t] * epsilon[i, t] * theta[i, t]
      
      # Likelihood
      Y[i, t] ~ dpois(mu[i, t])
    }
  }
  
  # Forward Filtering - inicialização em t = 1
  for (r in 1:n_regions) {
    a[r, 1] <- a0
    b[r, 1] <- b0
    
    g_b[r, 1] <- E[r, 1] * epsilon[r, 1] * exp(inprod(beta[1:n_covariates_theta], x[r, 1, 1:n_covariates_theta]))
    sum_g_b[r, 1] <- g_b[r, 1]
    sum_y[r, 1] <- count[r, 1]
    lambda[r,1] ~dgamma(shape = a[r,1], scale = 1/b[r,1])
    # lambda[r, 1] será atualizado fora do grafo (no sampler customizado)
  }
  
  # Forward Filtering para t >= 2
  for (t in 2:n_times) {
    for (r in 1:n_regions) {
      sum_y[r, t] <- sum_y[r, t - 1] + count[r, t]
      a[r, t] <- w * a[r, t - 1] + sum_y[r, t]
      
      g_b[r, t] <- E[r, t] * epsilon[r, t] * exp(inprod(beta[1:n_covariates_theta], x[r, t, 1:n_covariates_theta]))
      sum_g_b[r, t] <- sum_g_b[r, t - 1] + g_b[r, t]
      b[r, t] <- w * b[r, t - 1] + sum_g_b[r, t]
      
      lambda[r, t]~dgamma(shape = a[r,t], scale = 1/b[r,t])
    }
  }
})
#DEFININDO OS VALORES
modelConstants = list(K= K,
                      n_regions = dim(x)[1],
                      n_times = dim(x)[2],
                      n_covariates_theta = dim(x)[3], 
                      x = x,
                      hAI = hAI,
                      E = E,
                      sd_beta = 10,
                      a_unif = a_unif_ini,
                      b_unif = b_unif_ini,
                      a0 = a0,
                      b0 = b0,
                      mu_beta = beta_ini,
                      count =y_ini,
                      w=w
)
modelData = list(Y = round(y_ini))
modelInits1 = list(lambda = lambda_ini,
                  gamma = matrix(rep(gamma_ini,23),nrow =4,ncol =23),
                  beta = beta_ini)

modelInits2 = list(lambda =lambda_ini,gamma =matrix(rep(gamma_ini,23),nrow =4,ncol =23),beta =c(0,0,0) )
#Create the model
modelInits = c(modelInits1,modelInits2)
model = nimbleModel(code,constants = modelConstants,data = modelData,inits = modelInits)
#compilando o modelo
cModel = compileNimble(model)

#configurando os targets
modelConf = configureMCMC(model,monitors=c('beta','gamma',
                                           'mu',"lambda",'epsilon',
                                           'theta','a','b'), useConjugacy = TRUE,enableWAIC = F)

#removendo e adicionando o sampler
modelConf$removeSamplers("lambda")
modelConf$samplerConfs
modelConf$addSampler(target = "lambda",type = "dynamic_sampler",control = list(
  E = modelConstants$E,
  x = modelConstants$x,
  w = modelConstants$w,
  a0 = modelConstants$a0,
  b0 = modelConstants$b0,
  n_regions = modelConstants$n_regions,
  n_times = modelConstants$n_times
))
modelConf$samplerConfs

#construindo MCMC

modelMCMC = buildMCMC(modelConf)

cMCMC <- compileNimble(modelMCMC, project = model, resetFunctions = TRUE)

mcmc.out = runMCMC(cMCMC,
                      nchains = 2, niter = 100000,
                      summary = TRUE, WAIC = FALSE)

#Obter objetos mcmc separados por cadeia (para convergência)
samples_mcmc <- mcmc.out$samples  # lista de mcmc por cadeia

# Obter nomes de parâmetros
param_names <- varnames(samples_mcmc)

#traceplot dos betas
samples_beta1 = as.mcmc(samples_mcmc[,"beta[1]"])
traceplot(samples_beta1)
samples_beta2= as.mcmc(samples_mcmc[,"beta[2]"])
traceplot(samples_beta2)
samples_beta3 = as.mcmc(samples_mcmc[,"beta[3]"])
traceplot(samples_beta3)

samples_lambda1_1<-as.mcmc("lambda[1, 1]")
traceplot(samples_lambda1_1)
plot(samples_mcmc[,"lambda[1, 1]"],type =  "l")

# Filtrar nomes que NÃO sejam r[i] nem non.zero[i,t]
param_names_filtrados <- param_names[grepl("^lambda\\[[0-9]+, [0-9]+\\]$", param_names)]#&
#!grepl("^lambda\\[[0-9]+, [0-9]+\\]$", param_names)]


# Calcular a média a posteriori
media_posteriori <- colMeans(samples_mcmc)

#salvar cadeia

pdf(sprintf("C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/PesquisaMestrado/traceplots.pdf" ), width = 10, height = 8)
for (p in param_names_filtrados) {
  traceplot(samples_mcmc[, p], main = p)
  #abline(h = valores_reais[p], col = "red", lwd = 2)
}
dev.off()

media_posteriori
colMeans(mcmc.out$samples[["lambda[1, 1]"]])

mcmc_samples <- as.mcmc(mcmc.out$samples)

means_mcmc<-colMeans(mcmc_samples)

library(tidyverse)

# Convertendo o objeto MCMC para um data frame
mcmc_df <- as.data.frame(mcmc_samples)


# Transformando para formato long
mcmc_long <- pivot_longer(mcmc_df, everything(), names_to = "parameter", values_to = "value")

library(plotly)

# Criando um gráfico interativo de densidade
p <- ggplot(mcmc_long, aes(x = value, group = parameter)) +
  geom_density(aes(y = ..density..), fill = "blue", alpha = 0.5) +
  facet_wrap(~parameter, scales = "free") +
  labs(title = "Distribuição a Posteriori", x = "Valor", y = "Densidade")

# Convertendo para gráfico interativo
ggplotly(p)
