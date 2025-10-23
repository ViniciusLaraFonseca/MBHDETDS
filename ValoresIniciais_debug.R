# --- Script ValoresIniciais_Debug.R ---
source("_dataCaseStudy.R")
attach(data)
# Objetivo: Gerar valores iniciais E dados simulados controlados para depuração
#           do modelo MCMC dinâmico, seguindo a lógica de ValoresIniciais.R.

cat("--- Iniciando ValoresIniciais_Debug.R ---\n")

# --- PASSO 1: DEFINIR PARÂMETROS DA SIMULAÇÃO CONTROLADA ---

cat("--- PASSO 1: Definindo parâmetros 'verdadeiros' para a simulação ---\n")

set.seed(987) # Para reprodutibilidade da simulação

# Dimensões (use as dimensões reais do seu problema)
n_regions <- 75 
n_times   <- 23
p         <- 3 # Número de covariáveis beta
K         <- 4 # Número de clusters gamma

# Parâmetros "Verdadeiros" para Simulação
beta_true  <- c(-0.25, 0.5, -0.1) # Exemplo de valores verdadeiros
gamma_true <- c(0.05, 0.10, 0.10, 0.15)
w_true     <- 0.9
a0_true    <- 1 # Exemplo: um pouco mais informativo que 1
b0_true    <- 1 # Exemplo: rate

# --- PASSO 2: SIMULAR DADOS CONTROLADOS ---

cat("--- PASSO 2: Simulando dados controlados (Y, E, x) ---\n")

# Simular Covariáveis 'x'
x <- array(rnorm(n_regions * n_times * p), dim = c(n_regions, n_times, p))

# Simular e Normalizar Offset 'E'
E_raw <- matrix(runif(n_regions * n_times, 150, 250), nrow = n_regions) # Valores plausíveis
media_global_E_sim <- mean(E_raw)
E_debug <- E_raw #/ media_global_E_sim # Normalizado para média 1
cat("Offset 'E' simulado e normalizado (média 1).\n")

# Calcular Epsilon Fixo (baseado em gamma_true)
# (Assumindo que 'hAI' está carregado do _dataCaseStudy.r via attach(data))
if (!exists("hAI")) {
  stop("Erro: Objeto 'hAI' não encontrado. Carregue '_dataCaseStudy.r' e use attach(data).")
}
epsilon_ini <- matrix(NA, nrow=n_regions, ncol=1)
for(i in 1:n_regions) {
  epsilon_ini[i] <- 1 - sum(hAI[i, ] * gamma_true)
}


# Simular Trajetória "Verdadeira" de Lambda (usando a dinâmica do modelo)
# Inicialização de matrizes para o processo de simulação
lambda_at_fwd <- matrix(NA, ncol = n_times, nrow = n_regions)
y_at <- matrix(NA, nrow = n_regions, ncol = n_times + 1)
for(i in 1:n_regions){
  prod_val <- sum(x[i, 1, ] * beta_true)
  lambda0 <- rgamma(1,shape=1,rate=1)
  y_at[i, 1] <- E_debug[i,1]*lambda0*epsilon_ini[i]*exp(prod_val)
}
a0=1
b0=1
w=0.9
att_ini <- matrix(NA, nrow = n_regions, ncol = n_times)
btt_ini <- matrix(NA, nrow = n_regions, ncol = n_times)
at_ini  <- matrix(NA, nrow = n_regions, ncol = n_times + 1)
bt_ini  <- matrix(NA, nrow = n_regions, ncol = n_times + 1)

# Simulação Forward
for(i in 1:n_regions) {
  at_ini[i, 1] <- a0
  bt_ini[i, 1] <- b0
  for(t in 2:(n_times + 1)) {
    # Passos de filtragem (forward)
    att_ini[i, t-1] <- w * at_ini[i, t-1]
    btt_ini[i, t-1] <- w * bt_ini[i, t-1]
    at_ini[i, t]    <- att_ini[i, t-1] + y_at[i, t-1]
    
    prod_val <- sum(x[i, t-1, ] * beta_true)
    bt_ini[i, t]    <- btt_ini[i, t-1] + E_debug[i, t-1] * epsilon_ini[i] * exp(prod_val)
    
    # Amostra um valor de lambda a partir da distribuição filtrada (forward)
    lambda_at_fwd[i, t-1] <- rgamma(1, shape = at_ini[i, t], rate = bt_ini[i, t])
    
    # Simula um novo valor de Y para alimentar o próximo passo do filtro
    mu_at <- lambda_at_fwd[i, t-1] * exp(prod_val) * epsilon_ini[i] * E_debug[i, t-1]
    y_at[i, t] <- rpois(1, mu_at)
  }
}

# Simulação Backward (Suavização)
lambda_smooth <- matrix(NA, ncol = n_times, nrow = n_regions)
y_smooth <- matrix(NA, ncol = n_times, nrow = n_regions)

lambda_smooth[, n_times] <- lambda_at_fwd[, n_times] # O último ponto é o mesmo
y_smooth[, n_times] <- y_at[, n_times]

for(t in n_times:2) {
  # Amostra a "inovação" nu e calcula o lambda suavizado
  nu <- rgamma(n_regions, shape = (1-w) * at_ini[, t], rate = bt_ini[, t])
  lambda_smooth[, t-1] <- nu + w * lambda_smooth[, t]
  
  # Gera contagens Y consistentes com a trajetória suavizada de lambda
  prod_val_vec <- x[, t-1, ] %*% beta_true
  mu_smooth <- lambda_smooth[, t-1] * exp(prod_val_vec) * epsilon_ini * E_debug[, t-1]
  y_smooth[, t-1] <- rpois(n_regions, mu_smooth)
}

# --- PASSO 4: ATRIBUIÇÃO FINAL DOS VALORES INICIAIS ---

# Os objetos finais 'y_ini' e 'lambda_ini' contêm as trajetórias completas
# e temporalmente consistentes que serão usadas para iniciar o MCMC.
y_ini <- y_at[,2:24]
lambda_ini <- lambda_at_fwd

# Valores iniciais para os priors uniformes de gamma
a_unif_ini <- c(0.00, rep(0, K - 1))
b_unif_ini <- c(0.1, rep(0.99, K - 1))

beta_ini <- beta_true

# --- PASSO 4: ATRIBUIR VALORES FINAIS PARA USO NO MCMC ---

cat("--- PASSO 4: Atribuindo valores finais para o MCMC ---\n")

# ATENÇÃO: Sobrescrever os objetos globais com os valores simulados/estimados
Y <- y_ini      # Os dados observados para o modelo NIMBLE serão os simulados
E <- E_debug       # O offset para o modelo NIMBLE será o simulado/normalizado
# x já é x_debug, se você não carregou o Covariaveis.R original depois
# beta_ini já foi calculado
# gamma_ini já foi definido
# lambda_ini já foi calculado
a0 <- a0_true      # Hiperparâmetro para o modelo NIMBLE
b0 <- b0_true       # Hiperparâmetro para o modelo NIMBLE
w  <- w_true       # Fator de desconto para o modelo NIMBLE

cat("Valores iniciais e dados de depuração prontos.\n")
cat("Objetos criados/sobrescritos: Y, E, lambda_ini, beta_ini, gamma_ini, a0, b0, w\n")

# Limpeza de objetos intermédios (opcional)
# rm(list = ls()[grep("_debug|_ini$|_true$|_raw$|glm_fit_debug|lambda0_log_ini|mu_lambda0_ini|var_lambda0_ini", ls())])
# rm(media_global_E_sim, E_raw, epsilon_true, epsilon_matrix_true, lambda_true_sim, mu_true_sim)
# rm(lambda_at_fwd_ini, att_ini, btt_ini, at_sim, bt_sim, epsilon_ini_ffbs)