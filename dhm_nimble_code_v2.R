if(!require(dplyr)){ install.packages("dplyr"); require(dplyr)}   
if(!require(maps)){ install.packages("maps"); require(maps)}   # For adding a map to the plots of Brazil
if(!require(NbClust)){ install.packages("NbClust"); require(NbClust)}   
if(!require(ppclust)){ install.packages("ppclust"); require(ppclust)}   
if(!require(psych)){ install.packages("psych"); require(psych)}   
if(!require(nimble)){ install.packages("nimble"); require(nimble)} # For MCMC computation using NIMBLE.  
if(!require(spdep)){ install.packages("spdep"); require(spdep)}   # For computing the neighbourhood and adjancency objects.
if(!require(coda)){ install.packages("coda"); require(coda)}   
if(!require(viridis)){ install.packages("viridis"); require(viridis)}  # Colour pallette for plots. 
if(!require(ggplot2)){ install.packages("ggplot2"); require(ggplot2)}   


# Define the model
code <- nimbleCode({
  # Likelihood section
  for (t in 1:n_times) {  # Loop over time
    for (i in 1:n_regions) {  # Loop over regions
      # Calculation of epsilon (effect of covariates)
      epsilon[i, t] <- 1 - inprod(hAI[i, ], gamma[1:K, t])  # Inverse of the influence of gamma
      # Poisson likelihood
      Y[i, t] ~ dpois(mu[i, t])  # Observed counts Y for region i at time t
      
      # Mean of the Poisson distribution
      mu[i, t] <- lambda[i, t] * E[i, t] * epsilon[i, t] * theta[i, t]
      
      # Calculation of theta (exponential transformation of covariates)
      theta[i, t] <- exp(inprod(beta[1:n_covariates_theta], x[i, t, 1:n_covariates_theta]))
      
      # Terms for Poisson likelihood calculations
      #ppoi.term[i,t] <- exp(-mu[i, t] + Y[i, t] * log(mu[i, t] - lfactorial(Y[i, t])))
      #cpoi.term[i,t] <- 1 / ppoi.term[i,t]  # Cumulative Poisson term
    }
    
    # Prior for the first gamma parameter
    gamma[1, t] ~ dunif(min = a_unif[1], max = b_unif[1])
    
    # Priors for subsequent gamma parameters
    for (j in 2:K) {
      gamma[j, t] ~ dunif(min = a_unif[j] * (1 - sum(gamma[1:(j - 1), t])),
                          max = b_unif[j] * (1 - sum(gamma[1:(j - 1), t])))
    }
  }
  
  # Priors for beta coefficients
  for (l in 1:n_covariates_theta) {
    beta[l] ~ dnorm(mu_beta[l], sd_beta)  # Normal prior for each beta
  }
  # Estrutura foward lambda
  # Initialize time 1 properly
  for (r in 1:n_regions) {
    a[r, 1] <- a0
    b[r, 1] <- b0
    g_b[r, 1] <- E[r, 1] * epsilon[r, 1] * exp(inprod(beta[1:n_covariates_theta], x[r, 1, 1:n_covariates_theta]))
    sum_g_b[r, 1] <- g_b[r, 1]
    sum_y[r, 1] <- y_ini[r, 1]
    lambda[r, 1] ~ dgamma( shape = a[r, 1], scale = 1 / b[r, 1])
  }
  
  # Forward loop (t >= 2)
  for (t in 2:n_times) {
    for (r in 1:n_regions) {
      sum_y[r, t] <- sum_y[r, t-1] + y_ini[r, t]  # Recursive cumulative sum
      a[r, t] <- w * a[r, t-1] + sum_y[r, t]  # Use <- for assignment
      g_b[r, t] <- E[r, t] * epsilon[r, t] * exp(inprod(beta[1:n_covariates_theta] , x[r, t, 1:n_covariates_theta]))
      sum_g_b[r, t] <- sum_g_b[r, t-1] + g_b[r, t]  # Recursive cumulative
      b[r, t] <- w * b[r, t-1] + sum_g_b[r, t]  # Use <- for assignment
      lambda[r, t] ~ dgamma( shape = a[r, t], scale = 1 / b[r, t])  # rgamma not ~
    }
  }
})

#Define the constants 
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
                      w = w
);print(modelConstants)
#list with data
modelData = list(Y = round(y_ini));print(modelData)
#valor inicial tem que ta no espaço parametrico da priori
#list with initial values
modelInits = list(lambda = lambda_ini,
                  gamma = matrix(rep(gamma_ini,23),nrow =4,ncol =23),
                  beta = beta_ini);print(modelInits)
#Create the model

model = nimbleModel(code,constants = modelConstants,data = modelData,inits = modelInits)

model$defaultModelValues
#checking the model
print(model)
print(model$lambda)
print(model$gamma)
print(model$beta)
print(model$epsilon)
model$calculate()
cat("Dimensions of x:", dim(x), "\n")
cat("Dimensions of Y:", dim(y_ini), "\n")
#compilando o modelo
cModel = compileNimble(model)

#configurando os targets
modelConf = configureMCMC(model,monitors=c('Y','beta','gamma',
                                           'mu',"lambda",'epsilon',
                                           'theta','a','b'), useConjugacy = TRUE,enableWAIC = TRUE)

#removendo e adicionando o sampler
modelConf$removeSamplers("lambda")
modelConf$addSampler(target = "lambda",type = "dynamic_sampler")
modelConf$samplerConfs


#construindo MCMC

modelMCMC = buildMCMC(modelConf)

cMCMC <- compileNimble(modelMCMC, project = model, resetFunctions = TRUE)

ls()
#MCMC

mcmc.out = nimbleMCMC(cMCMC,inits = modelInits,
                      nchains = 1, niter = 10000,
                      summary = TRUE, WAIC = FALSE)

mcmc.out$samples
gelman.diag(mcmc.out$samples[,c('beta[1]','beta[2]','beta[3]')])
aux = mcmc.out$samples$chain1
plot(aux[,1])  
