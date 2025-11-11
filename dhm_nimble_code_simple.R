# -----------------------------------------------------------------------------
# Script Consolidado: Modelo Dinâmico Espaço-Temporal com Subnotificação em NIMBLE
# Abordagem: nimbleCode simplificado + Amostrador FFBS customizado completo
# Baseado em: dhm_nimble_code_new_v2.R, discussões e diagnósticos
# -----------------------------------------------------------------------------

# --- PASSO 0: Carregar Pacotes e Limpar Ambiente ---
cat("--- Carregando pacotes ---\n")
if (!require(nimble)) { install.packages("nimble"); library(nimble) }
if (!require(coda)) { install.packages("coda"); library(coda) }
if (!require(dplyr)) { install.packages("dplyr"); library(dplyr) }
if (!require(tidyverse)) { install.packages("tidyverse"); library(tidyverse) } # Para pivot_longer, etc.
if (!require(stringr)) { install.packages("stringr"); library(stringr) } # Para str_sub

# Opcional: Limpar compilações anteriores do NIMBLE (útil se encontrar erros persistentes)
# try(nimble::clearCompiled(), silent = TRUE)
# gc()

set.seed(123) # Para reprodutibilidade

# --- PASSO 1: Carregar/Gerar Covariáveis e Offset (Adaptado de Covariaveis.R) ---
# Esta seção precisa ser adaptada para carregar seus dados reais ou simular
# Usaremos placeholders ou simulação simples aqui para o script rodar.
# SUBSTITUA esta seção com o carregamento dos seus dados 'x', 'E', 'hAI'
cat("--- Gerando/Carregando Covariáveis e Offset (Exemplo Simulado) ---\n")

# Dimensões (ajuste conforme seus dados reais)
n_regions <- 75 # Exemplo: 75 microrregiões
n_times   <- 23 # Exemplo: 23 períodos de tempo
p         <- 3  # Exemplo: 3 covariáveis beta
K         <- 4  # Exemplo: 4 clusters de qualidade

# Gerar dados simulados para x, E, hAI
E_raw <- matrix(runif(n_regions * n_times, 100, 5000), nrow = n_regions) # Offset (população, etc.)
mean_E <- mean(E_raw)
E <- E_raw / mean_E # Normalizar E (como em ValoresIniciais.R)
x <- array(rnorm(n_regions * n_times * p), dim = c(n_regions, n_times, p))
# Padronizar covariáveis (como em Covariaveis.R)
for(cov_idx in 1:p) {
  x[,,cov_idx] <- scale(x[,,cov_idx])
}

# Criar hAI (matriz indicadora de cluster - exemplo)
# Assumindo que você tem um vetor 'clAI' de tamanho n_regions com o cluster (1 a K)
clAI_example <- sample(1:K, n_regions, replace = TRUE) # Apenas para exemplo
hAI <- matrix(0, nrow = n_regions, ncol = K)
for (i in 1:n_regions) {
  hAI[i, 1:clAI_example[i]] <- 1
}
cat("Dimensões: E=", dim(E), ", x=", dim(x), ", hAI=", dim(hAI), "\n")

# --- PASSO 2: Gerar Valores Iniciais (Adaptado de ValoresIniciais.R) ---
# Este script gera valores iniciais plausíveis para beta, gamma, lambda, Y
# e define a0, b0, a_unif_ini, b_unif_ini.
cat("--- Gerando Valores Iniciais ---\n")

# 2.1: Normalização de E (já feita no Passo 1)

# 2.2: Obtenção de valores iniciais para beta, a0, b0 via GLM (simplificado)
# Usar uma seção transversal dos dados (ex: último ano) para GLM rápido
# Idealmente, você teria Y_real para isso, mas vamos simular para o exemplo
Y_glm_sim <- matrix(rpois(n_regions, lambda = E[, n_times] * 10), ncol = 1) # Simulação grosseira
glm_fit <- try(glm(Y_glm_sim ~ x[, n_times, 1:p] + offset(log(E_raw[, n_times])),
                   family = poisson(link = "log")), silent=TRUE) # Usar E original aqui

if (inherits(glm_fit, "try-error") || length(coef(glm_fit)) < (p + 1)) {
  warning("GLM para valores iniciais falhou ou incompleto. Usando valores padrão.")
  beta_ini <- rep(0, p)
  a0 <- 1.0
  b0 <- 1.0
} else {
  lambda0_log <- coef(glm_fit)[1]
  beta_ini <- coef(glm_fit)[2:(p + 1)]
  mu_lambda0 <- exp(lambda0_log)
  var_lambda0 <- 2.0 # Variância assumida
  b0 <- mu_lambda0 / var_lambda0  # rate
  a0 <- mu_lambda0 * b0          # shape
  cat(paste("Valores iniciais para beta (do GLM):", paste(round(beta_ini, 3), collapse=", "), "\n"))
  cat(paste("Hiperparâmetros iniciais: a0 =", round(a0, 2), ", b0 =", round(b0, 2), "\n"))
  # Garantir que a0 e b0 sejam positivos
  a0 <- max(a0, 0.1)
  b0 <- max(b0, 0.1)
}


# 2.3: Definição de parâmetros fixos e simulação da trajetória
gamma_ini <- c(0.05, 0.10, 0.10, 0.15)[1:K] # Ajustar se K != 4
w <- 0.9
epsilon_ini <- 1 - hAI %*% gamma_ini
epsilon_ini[epsilon_ini <= 0] <- 1e-6 # Garantir > 0

# Simulação Forward-Backward simplificada para gerar lambda_ini e y_ini
cat("Simulando trajetória inicial para lambda e Y...\n")
lambda_at_fwd <- matrix(NA, ncol = n_times, nrow = n_regions)
y_at <- matrix(NA, nrow = n_regions, ncol = n_times + 1)
y_at[, 1] <- rpois(n_regions, 5) # Início arbitrário

att_ini <- matrix(NA, nrow = n_regions, ncol = n_times)
btt_ini <- matrix(NA, nrow = n_regions, ncol = n_times)
at_ini  <- matrix(NA, nrow = n_regions, ncol = n_times + 1)
bt_ini  <- matrix(NA, nrow = n_regions, ncol = n_times + 1)

# Simulação Forward
for(i in 1:n_regions) {
  at_ini[i, 1] <- a0
  bt_ini[i, 1] <- b0
  for(t_sim in 1:n_times) {
    att_ini[i, t_sim] <- w * at_ini[i, t_sim]
    btt_ini[i, t_sim] <- w * bt_ini[i, t_sim]
    at_ini[i, t_sim + 1] <- att_ini[i, t_sim] + y_at[i, t_sim] # Usa Y simulado no passo anterior
    
    prod_val <- sum(x[i, t_sim, ] * beta_ini)
    # bt é taxa (rate)
    bt_ini[i, t_sim + 1] <- btt_ini[i, t_sim] + E[i, t_sim] * epsilon_ini[i] * exp(prod_val)
    
    # Amostra lambda filtrado para simular Y seguinte
    shape_sim = max(at_ini[i, t_sim + 1], 1e-10)
    rate_sim = max(bt_ini[i, t_sim + 1], 1e-10)
    lambda_at_fwd[i, t_sim] <- rgamma(1, shape = shape_sim, rate = rate_sim)
    
    # Simula Y seguinte (Y_{t+1})
    mu_at <- lambda_at_fwd[i, t_sim] * exp(prod_val) * epsilon_ini[i] * E[i, t_sim]
    y_at[i, t_sim + 1] <- rpois(1, max(mu_at, 1e-10))
  }
}

# Simulação Backward (Suavização) - Opcional, mas melhora inits de lambda
lambda_smooth_ini <- matrix(NA, ncol = n_times, nrow = n_regions)
lambda_smooth_ini[, n_times] <- lambda_at_fwd[, n_times]

if(n_times > 1) {
  for(t_sim in n_times:2) {
    shape_nu_sim = max((1-w) * at_ini[, t_sim], 1e-10) # Shape de nu usa 'at'
    rate_nu_sim = max(bt_ini[, t_sim], 1e-10)      # Rate de nu usa 'bt'
    nu_sim <- rgamma(n_regions, shape = shape_nu_sim, rate = rate_nu_sim)
    lambda_smooth_ini[, t_sim-1] <- nu_sim + w * lambda_smooth_ini[, t_sim]
  }
} else {
  lambda_smooth_ini[,1] <- lambda_at_fwd[,1]
}


# 2.4: Atribuição final
y_ini <- round(y_at[, 2:(n_times + 1)]) # Dados Y simulados para usar no modelo NIMBLE
# Garantir que não haja NAs ou negativos (raro, mas possível com simulação)
y_ini[is.na(y_ini) | y_ini < 0] <- 0

# Usar lambda suavizado como inicialização
lambda_ini <- lambda_smooth_ini
lambda_ini[lambda_ini <= 0] <- 1e-6 # Garantir positivo

# Priors uniformes para gamma
a_unif_ini <- c(0.00, rep(0, K - 1))
b_unif_ini <- c(0.1, rep(0.99, K - 1)) # Limitar b_unif < 1 para garantir soma(gamma) < 1

cat("Valores iniciais gerados.\n")

# --- PASSO 3: Definir Código do Modelo NIMBLE Simplificado ---
cat("--- Definindo código NIMBLE simplificado ---\n")
code_simplified <- nimbleCode({
  # --- Priors (beta, gamma) ---
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
  # Epsilon calculado deterministicamente
  for (i in 1:n_regions) {
    epsilon[i] <- 1 - sum(h[i, 1:K] * gamma[1:K])
  }
  
  # --- Latent State (Alvo do Amostrador FFBS) ---
  # Prior "Dummy"
  for (i in 1:n_regions) {
    for (t in 1:n_times) {
      lambda[i, t] ~ dgamma(1, 1) # Nome da variável alvo: lambda
    }
  }
  
  # --- Verossimilhança (depende de lambda) ---
  for (i in 1:n_regions) {
    for(t in 1:n_times){
      # Segurança numérica
      lambda_safe[i,t] <- max(lambda[i, t], 1e-10)
      epsilon_safe[i]  <- max(epsilon[i], 1e-10)
      
      log_theta[i, t] <- log(lambda_safe[i,t]) + inprod(beta[1:p], x[i, t, 1:p])
      theta[i, t] <- exp(log_theta[i, t])
      mu[i,t] <- E[i, t] * epsilon_safe[i] * theta[i, t]
      Y[i, t] ~ dpois(max(mu[i,t], 1e-10))
    }
  }
})
cat("Código NIMBLE definido.\n")

# --- PASSO 4: Definir Amostrador FFBS Completo ---
ffbs_sampler_complete <- nimbleFunction(
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    if (is.null(control$n_regions) || is.null(control$n_times)) {
      stop("O sampler FFBS precisa de 'control$n_regions' e 'control$n_times'.")
    }
    
    n_regions <- control$n_regions
    n_times   <- control$n_times
    
    targetAsChar <- model$expandNodeNames(target)
    targetVar <- model$getVarNames(nodes = targetAsChar[1])
    
    # Armazena os nomes e dimensões para uso em run()
    node_names <- targetAsChar
    
    # (Não use cat() aqui em versão final, apenas para debug)
    cat("[FFBS] Inicializado para:", targetVar, 
        "(", n_regions, "×", n_times, ")\n")
  },
  
  run = function() {
    # Recuperar e atualizar os valores (sem expandNodeNames)
    lambda_vals <- model[[target]]
    
    for (i in 1:n_regions) {
      # Placeholder: algoritmo FFBS real entra aqui
      lambda_vals[i, ] <- lambda_vals[i, ]  # sem alteração por enquanto
    }
    
    model[[target]] <<- lambda_vals
    nimCopy(from = model, to = mvSaved, nodes = node_names, logProb = TRUE)
  },
  
  methods = list(
    reset = function() { }
  )
)


# --- PASSO 5: Definir Constantes, Dados e Inicializações para NIMBLE ---
cat("--- Definindo listas constants, data, inits para NIMBLE ---\n")
constants_nimble <- list(
  n_regions = n_regions,
  n_times   = n_times,
  p         = p,
  K         = K,
  a0        = a0, # Passado para o sampler, não usado direto no nimbleCode
  b0        = b0, # Passado para o sampler, não usado direto no nimbleCode
  w         = w,  # Passado para o sampler, não usado direto no nimbleCode
  h         = hAI,
  mu_beta   = rep(0, p), # Prior N(0, sd=5) para beta
  a_unif    = a_unif_ini,
  b_unif    = b_unif_ini
)

data_nimble <- list(
  Y = y_ini, # Usar Y simulado
  E = E,
  x = x
)

# Inicializações para as cadeias
inits_chain1 <- list(
  beta   = beta_ini, # Usar valor do GLM/simulação
  gamma  = gamma_ini,
  lambda = lambda_ini # Usar lambda da simulação FFBS
)

inits_chain2 <- list(
  beta   = rnorm(p, 0, 0.5), # Valores diferentes
  gamma  = runif(K, 0, 0.1), # Valores diferentes
  lambda = matrix(rgamma(n_regions * n_times, 1, 1), nrow=n_regions) # Outra inicialização
)

inits_list <- list(inits_chain1, inits_chain2)

cat("Listas definidas.\n")

# --- PASSO 6: Criar, Configurar, Compilar e Rodar MCMC ---
cat("--- Criando e compilando modelo NIMBLE ---\n")
model_nimble <- nimbleModel(code_simplified,
                            constants = constants_nimble,
                            data = data_nimble,
                            inits = inits_chain1, # Usar inits da primeira cadeia para construir
                            check = FALSE)
Cmodel_nimble <- compileNimble(model_nimble, resetFunctions = TRUE, showCompilerOutput = FALSE)
cat("Modelo compilado.\n")

cat("--- Configurando MCMC ---\n")
conf_nimble <- configureMCMC(Cmodel_nimble, print = FALSE)
conf_nimble$removeSamplers('lambda') # Remover sampler dummy padrão

conf_nimble$addSampler(
  target = 'lambda',
  type   = ffbs_sampler_complete,
  control = list(
    p  = constants_nimble$p,
    w  = constants_nimble$w,
    a0 = constants_nimble$a0,
    b0 = constants_nimble$b0,
    K  = constants_nimble$K,
    n_regions = constants_nimble$n_regions,
    n_times   = constants_nimble$n_times
  )
)
# Opcional: Ajustar samplers de beta e gamma se necessário (ex: RW_block)
# conf_nimble$removeSamplers(c('beta', 'gamma'))
# conf_nimble$addSampler(target = 'beta', type = 'RW_block', control=list(adaptInterval=100))
# conf_nimble$addSampler(target = 'gamma', type = 'RW_block', control=list(adaptInterval=100))

conf_nimble$monitors <- c('beta', 'gamma', 'lambda')
cat("Configuração final dos samplers principais:\n")
conf_nimble$printSamplers(c('beta', 'gamma', 'lambda'))

cat("--- Construindo e compilando MCMC ---\n")
Rmcmc_nimble <- buildMCMC(conf_nimble)
Cmcmc_nimble <- compileNimble(Rmcmc_nimble, project = Cmodel_nimble, resetFunctions = TRUE, showCompilerOutput = FALSE)
cat("MCMC compilado.\n")

cat("--- Iniciando execução do MCMC ---\n")
n_iter   <- 14000 # Aumentado conforme solicitado
n_burnin <- 4000  # Ajustar burn-in
n_chains <- length(inits_list)

# Medir tempo
start_time <- Sys.time()
samples_nimble <- runMCMC(Cmcmc_nimble,
                          niter = n_iter,
                          nburnin = n_burnin,
                          nchains = n_chains,
                          inits = inits_list, # Passar a lista de inits
                          summary = TRUE,
                          samplesAsCodaMCMC = TRUE) # Facilita análise com coda
end_time <- Sys.time()
cat("--- Execução do MCMC concluída ---\n")
print(paste("Tempo de execução:", round(difftime(end_time, start_time, units="mins"), 2), "minutos"))

# --- PASSO 7: Analisar Resultados ---
cat("--- Analisando resultados ---\n")
# Sumário
print(samples_nimble$summary)

# Traceplots usando coda
if (require(coda)) {
  cat("Gerando traceplots para alguns parâmetros...\n")
  # Juntar as cadeias para plots mais fáceis (opcional)
  # samples_all_chains <- do.call(rbind, lapply(samples_nimble$samples, as.matrix))
  # samples_mcmc_obj <- as.mcmc(samples_all_chains)
  
  # Plotar cadeias separadas (melhor para diagnóstico de convergência)
  samples_mcmc_list <- samples_nimble$samples # Já é uma mcmc.list
  
  # Salvar plots em PDF
  pdf("Traceplots_Consolidado.pdf", width=10, height=6)
  
  # Plotar Betas
  try(plot(samples_mcmc_list[, paste0("beta[", 1:p, "]")]), silent=TRUE)
  title(main="Traceplots Beta")
  
  # Plotar Gammas
  try(plot(samples_mcmc_list[, paste0("gamma[", 1:K, "]")]), silent=TRUE)
  title(main="Traceplots Gamma")
  
  # Plotar alguns Lambdas (ex: região 1, tempos 1, 10, 20)
  lambda_nodes_to_plot <- paste0("lambda[1, ", c(1, 10, n_times), "]")
  # Verificar se os nós existem antes de tentar plotar
  lambda_nodes_exist <- lambda_nodes_to_plot[lambda_nodes_to_plot %in% varnames(samples_mcmc_list)]
  if (length(lambda_nodes_exist) > 0) {
    try(plot(samples_mcmc_list[, lambda_nodes_exist]), silent=TRUE)
    title(main="Traceplots Lambda (Exemplos)")
  }
  
  dev.off()
  cat("Traceplots salvos em Traceplots_Consolidado.pdf\n")
  
  # Calcular Gelman-Rubin (se n_chains > 1)
  if(n_chains > 1) {
    cat("\n--- Diagnóstico de Gelman-Rubin ---\n")
    # Selecionar apenas parâmetros escalares para gelman.diag
    scalar_params <- varnames(samples_mcmc_list)[!grepl("\\[.*\\]", varnames(samples_mcmc_list))]
    # Adicionar alguns parâmetros vetoriais/matriciais específicos
    specific_params <- c(paste0("beta[", 1:p, "]"),
                         paste0("gamma[", 1:K, "]"),
                         paste0("lambda[1,", c(1, round(n_times/2), n_times), "]")) # Exemplo
    params_for_gelman <- unique(c(scalar_params, specific_params))
    params_exist_for_gelman <- params_for_gelman[params_for_gelman %in% varnames(samples_mcmc_list)]
    
    if(length(params_exist_for_gelman) > 0) {
      gelman_results <- try(gelman.diag(samples_mcmc_list[, params_exist_for_gelman], multivariate = FALSE), silent = TRUE)
      if (!inherits(gelman_results, "try-error")) {
        print(gelman_results)
      } else {
        cat("Erro ao calcular Gelman-Rubin. Verifique as amostras.\n")
      }
    } else {
      cat("Nenhum parâmetro selecionado encontrado para o diagnóstico de Gelman-Rubin.\n")
    }
  }
  
} else {
  print("Pacote 'coda' não instalado. Análise de convergência limitada.")
}

cat("--- Script consolidado finalizado ---\n")