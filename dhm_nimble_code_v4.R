library(stringr)

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
a0 = 2.15  # Parâmetro a inicial
b0 = 0.464  # Parâmetro b inicial

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
lambda_ini = matrix(rgamma(1725, shape = a0, rate = b0), ncol = 23, nrow = 75)

#estimando o theta com os betas iniciais e o x
theta_ini = lambda_ini * exp(beta_ini[1] * x[,, 1] )

# Estimando Y com base nos parâmetros lambda, E e epsilon, aplicando a função exponencial nos coeficientes beta
y_ini =  E * epsilon_ini * theta_ini

#estimando o theta com os betas iniciais e o x
theta_ini = lambda_ini*exp(beta_ini * x[,, 1] )

# Definindo um parâmetro w para ser utilizado em cálculos posteriores
w=0.7

# código nimble

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

#DEFININDO OS VALORES
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
modelData = list(Y = round(y_ini),
                 x = x[,,1],
                 hAI = hAI,
                 E = E)
modelInits = list(lambda = lambda_ini,
                   gamma = gamma_ini,
                   beta = beta_ini)


#Create the model
#modelInits = c(modelInits1,modelInits2)
model = nimbleModel(code,constants = modelConstants,data = modelData,inits = modelInits)
#compilando o modelo
cModel = compileNimble(model)

#configurando os targets
modelConf = configureMCMC(model,monitors=c('beta','gamma',
                                           'mu',"lambda",'epsilon',
                                           'theta'), useConjugacy = TRUE,enableWAIC = F)

#removendo e adicionando o sampler
modelConf$removeSamplers("lambda")
modelConf$samplerConfs
modelConf$addSampler(target = "lambda",type = "dynamic_sampler",control = list(
  E = modelData$E,
  x = modelData$x,
  w = modelConstants$w,
  a0 = modelConstants$a0,
  b0 = modelConstants$b0,
  n_regions = modelConstants$n_regions,
  n_times = modelConstants$n_times
))
modelConf$removeSampler("beta")
modelConf$addSampler("beta", type = "RW_block",
                    control = list(scale = 0.1, adaptScale = TRUE))
#modelConf$addSampler(target = "beta",type = "slice")
modelConf$samplerConfs
#construindo MCMC

modelMCMC = buildMCMC(modelConf)

cMCMC <- compileNimble(modelMCMC, project = model, resetFunctions = TRUE)

mcmc.out = runMCMC(cMCMC,
                   nchains = 2, niter = 50000,
                   summary = TRUE, WAIC = FALSE)
samples <-as.mcmc(mcmc.out$samples)
chain1 = as.mcmc(mcmc.out$samples$chain1)
chain2 = as.mcmc(mcmc.out$samples$chain2)

mcmc_sample=mcmc.list(chain1,chain2)

#traceplot betas
traceplot(chain1[,"beta"])
traceplot(chain2[,"beta"])

#traceplot theta
traceplot(chain1[,"theta[1, 1]"])
traceplot(chain2[,"theta[1, 1]"])

#traceplot lambda
traceplot(chain1[,"lambda[1, 1]"],ylim = c(0,30))
traceplot(chain2[,"lambda[1, 10]"])
traceplot(chain2[,"lambda[1, 23]"])
#traceplot theta
traceplot(chain1[,"gamma[1]"])
traceplot(chain1[,"gamma[2]"])
traceplot(chain1[,"gamma[3]"])
traceplot(chain1[,"gamma[4]"])

#traceplot  mu
traceplot(chain1[,'mu[1, 1]'])
traceplot(chain1[,'mu[1, 23]'])
#manipulção dos dados para analise gráfica das médias a posteriori

dfchain1 = as.data.frame(chain1)
dfchain1 = colMeans(dfchain1)
dfchain1= as.data.frame(dfchain1)
dfchain1_long <- dfchain1 %>%
  rownames_to_column(var = "var") %>%
  mutate(
    base_var = str_extract(var, "^(\\w+)"),
    region = as.integer(str_extract(var, "(?<=\\[)(\\d+)(?=,)")),
    time = as.integer(str_extract(var, "(?<=, )(\\d+)(?=\\])"))
  )
varnames(mcmc.out$samples)

dfchain1_long%>%
  filter(base_var == "lambda")%>%
  filter(region == 1)%>%
  ggplot(aes(x = time,y=dfchain1 ))+geom_line()

dfchain1_long%>%
  filter(base_var == "lambda")%>%
  filter(region %in% c(1,30,75))%>%
  ggplot(aes(x = time,y=dfchain1, color = as.factor(region)))+geom_line()

dfchain2 = as.data.frame(chain2)
dfchain2 = colMeans(dfchain2)
dfchain2 = as.data.frame(dfchain2)
dfchain2_long <- dfchain2 %>%
  rownames_to_column(var = "var") %>%
  mutate(
    base_var = str_extract(var, "^(\\w+)"),
    region = as.integer(str_extract(var, "(?<=\\[)(\\d+)(?=,)")),
    time = as.integer(str_extract(var, "(?<=, )(\\d+)(?=\\])"))
  )
varnames(mcmc.out$samples)

dfchain2_long%>%
  filter(base_var == "lambda")%>%
  filter(region == 1)%>%
  ggplot(aes(x = time,y=dfchain1 ))+geom_line()


dfchain2_long%>%
  filter(base_var == "lambda")%>%
  filter(region %in% c(1,30,75))%>%
  ggplot(aes(x = time,y=dfchain2, color = as.factor(region)))+geom_line()
