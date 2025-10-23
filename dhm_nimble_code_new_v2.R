# --- PASSO 0: CARREGAR PACOTES E GARANTIR REPRODUTIBILIDADE ---
if (!require(nimble)) install.packages("nimble")
if (!require(coda)) install.packages("coda")
if (!require(ggplot2)) install.packages("ggplot2") # Para os gráficos
library(nimble)
library(coda)
library(ggplot2)

set.seed(456) # Para reprodutibilidade

cat("--- Início do Script de Depuração com Traceplots ---\n\n")

# --- PASSO 1: CARREGAR DADOS, COVARIÁVEIS E VALORES INICIAIS ---
cat("--- PASSO 1: Carregando dados e valores iniciais ---\n")
# (Certifique-se que a linha com 'require(fda)' em Covariaveis.R está comentada)
source("_dataCaseStudy.r")
attach(data)
source("Covariaveis.R")
source("ValoresIniciais_debug.R")
# Guardar os valores "verdadeiros" da simulação para plotar
lambda_verdadeiro <- lambda_ini
beta_verdadeiro <- beta_ini
gamma_verdadeiro <- gamma_ini

# --- PASSO 2: DEFINIR O MODELO ESTILO v8 ---
cat("--- PASSO 2: Definindo o modelo estilo v8 ---\n")

code_v8_style <- nimbleCode({
  # Priors (beta e gamma são amostrados)
  for (j in 1:p) { beta[j] ~ dnorm(mu_beta[j], sd = 5) }
  gamma[1] ~ dunif(a_unif[1], b_unif[1])
  for (j in 2:K) {
    gamma[j] ~ dunif(
      min = a_unif[j] * (1 - sum(gamma[1:(j - 1)])),
      max = b_unif[j] * (1 - sum(gamma[1:(j - 1)]))
    )
  }
  for (i in 1:n_regions) { epsilon[i] <- 1 - sum(h[i, 1:K] * gamma[1:K]) }
  
  # Lógica do Filtro Forward
  for (i in 1:n_regions) {
    at[i, 1] <- a0
    bt[i, 1] <- b0
    
    for (t in 2:(n_times + 1)) {
      att[i, t-1] <- w * at[i, t-1]
      btt[i, t-1] <- w * bt[i, t-1]
      at[i, t] <- att[i, t-1] + Y[i, t-1]
      bt[i, t] <- btt[i, t-1] + E[i, t-1] * epsilon[i] * exp(inprod(beta[1:p], x[i, t-1, 1:p]))
      #DISTRIBUIÇÃO TERMO DE INOVAÇÃO (IMPORTANTE PRO PASSO BACKWARD)
      nu[i,t-1] ~dgamma((1-w)at[i,t-1],bt[i,t-1])
    }
    #PREDITIVA PRIORI DO LAMBDA 
    for (t in 1:n_times) {
      lambda[i, t] ~ dgamma(att[i, t], rate = btt[i, t])
    }
    lambda_suv[i,n_times] <- lambda[i,n_times]
    #Lógica do backward
    for(t in 1:(n_times-1)) {
      lambda_suv[i, n_times-t] <- nu[i,n_times-t] + w * lambda[i,n_times-t+1]
    }
    #AGORA COM TODAS AS PEÇAS, VAMOS PARA A VEROSSIMILHANÇA
    for (t in 1:n_times) {
      mu[i,t] <- E[i, t] * epsilon[i] * lambda_suv[i, t] * exp(inprod(beta[1:p], x[i, t, 1:p]))
      Y[i, t] ~ dpois(mu[i,t])
    }
  }
})

# --- PASSO 3: CONFIGURAR E EXECUTAR O MCMC (COMPILADO) ---
cat("\n--- PASSO 3: Configurando e executando o MCMC ---\n")

n_regions <- dim(E)[1]; n_times <- dim(E)[2]; p_params <- dim(x)[3]; K_clusters <- 4
constants_full <- list(n_regions = n_regions, n_times = n_times, p = p_params, K = K_clusters, h = hAI, 
                       mu_beta = beta_ini, a_unif = a_unif_ini, b_unif = b_unif_ini, a0=a0, b0=b0, w=w)
data_full <- list(Y = Y, E = E, x = x)

# Valores iniciais para duas cadeias
inits1 <- list(lambda = lambda_ini, beta = beta_ini, gamma = gamma_ini)
inits2 <- list(lambda = matrix(1, n_regions, n_times), beta = rnorm(p_params, 0, 0.5), gamma = gamma_ini / 2)
inits_list <- list(inits1, inits2)

# Construir e compilar o modelo
model_v8 <- nimbleModel(code_v8_style, constants = constants_full, data = data_full, inits = inits1, check = FALSE)
Cmodel_v8 <- compileNimble(model_v8)

# Configurar MCMC (NIMBLE usará seus amostradores padrão, o que é o objetivo do teste)
conf_v8 <- configureMCMC(Cmodel_v8, monitors = c("lambda", "beta", "gamma","lambda_suv","nu"))
Rmcmc_v8 <- buildMCMC(conf_v8)
Cmcmc_v8 <- compileNimble(Rmcmc_v8, project = Cmodel_v8)

# Executar o MCMC com 2 cadeias
cat("\n--- Executando MCMC com 2 cadeias (pode demorar alguns minutos) ---\n")
samples <- runMCMC(Cmcmc_v8, niter = 5000, nburnin = 1000, nchains = 2, inits = inits_list, samplesAsCodaMCMC = TRUE)

# --- PASSO 4: GERAR E SALVAR OS TRACEPLOTS ---
cat("\n--- PASSO 4: Gerando e salvando os traceplots ---\n")

# Converter para um formato mais fácil de plotar com ggplot
samples_df <- do.call(rbind, lapply(1:length(samples), function(i) {
  df <- as.data.frame(samples[[i]])
  df$chain <- as.factor(i)
  df$iter <- 1:nrow(df)
  return(df)
}))

# Função para criar e salvar traceplots
save_traceplot <- function(df, param_name, true_value) {
  p <- ggplot(df, aes_string(x = "iter", y = paste0("`", param_name, "`"), color = "chain")) +
    geom_line(alpha = 0.7) +
    geom_hline(yintercept = true_value, color = "blue", linetype = "dashed", size = 1) +
    labs(title = paste("Traceplot for", param_name), 
         subtitle = paste("Blue dashed line = True simulated value (", round(true_value, 3), ")"),
         x = "Iteration", y = "Value") +
    theme_minimal() +
    scale_color_manual(values = c("black", "red"))
  
  filename <- paste0("traceplot_v8_", gsub("\\[|, |\\]", "_", param_name), ".png")
  ggsave(filename, p, width = 8, height = 4)
  cat(paste("Salvo:", filename, "\n"))
}

# Gerar plots para beta
for(i in 1:p_params) {
  param <- paste0("beta[", i, "]")
  save_traceplot(samples_df, param, beta_verdadeiro[i])
}

# Gerar plots para gamma
for(i in 1:K_clusters) {
  param <- paste0("gamma[", i, "]")
  save_traceplot(samples_df, param, gamma_verdadeiro[i])
}

# Gerar plots para lambdas selecionados
selected_lambdas <- c("lambda[1, 1]", "lambda[1, 10]", "lambda[1, 23]")
for(lam_name in selected_lambdas) {
  # Extrair os índices para obter o valor verdadeiro
  indices <- as.numeric(unlist(regmatches(lam_name, gregexpr("[0-9]+", lam_name))))
  true_val <- lambda_verdadeiro[indices[1], indices[2]]
  save_traceplot(samples_df, lam_name, true_val)
}

cat("\n--- Análise concluída. Verifique os ficheiros .png no seu diretório. ---\n")